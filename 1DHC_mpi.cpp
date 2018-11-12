/*
--1D heat conduction 
--FT using ULMF
 */
# include <cmath>
# include <cstdlib>
# include <ctime>
# include <fstream>
# include <iostream>
# include <mpi.h>
# include <mpi-ext.h>
# include <signal.h>
# include <setjmp.h>

using namespace std;

int main ( int argc, char *argv[] );
void htc ( int rank, int np );
static void verbose_errhandler(MPI_Comm* pcomm, int* perr, ...);
static int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm);
static int app_buddy_ckpt(MPI_Comm comm);
static int app_reload_ckpt(MPI_Comm comm);
static int app_needs_repair(void);
static int iteration = 0;
static int ckpt_iteration = -1;
static const int ckpt_tag = 42;
static jmp_buf restart;
#define lbuddy(r) ((r+np-1)%np)
#define rbuddy(r) ((r+np+1)%np)

static MPI_Comm world = MPI_COMM_NULL;

//****************************************************************************

int main ( int argc, char *argv[] )

//****************************************************************************
{
  int rank;
  int size, np;
  int victim;
  MPI_Errhandler errh;
  MPI_COMM parent;

//
//  Initialize MPI.
//

  MPI_Init ( &argc, &argv );
  MPI_Comm_create_errhandler(&errhandler_respawn, &errh);

    /* 1st time or repair? */
  MPI_Comm_get_parent( &parent );
  if( MPI_COMM_NULL == parent )
    {
        /* First run: Let's create an initial world,
         * a copy of MPI_COMM_WORLD */
        MPI_Comm_dup( MPI_COMM_WORLD, &world );
        MPI_Comm_size( world, &np );
        MPI_Comm_rank( world, &rank );
    }

  else
    {
        /* repair, lets get the repaired world */
        app_needs_repair();
    }
  /* We set an errhandler on world, so that a failure is not fatal anymore. */
  MPI_Comm_set_errhandler( world, errh );
    
  htc ( rank, np );
  
//
//  Terminate MPI.
//
  MPI_Finalize ( );

  return 0;
}
//****************************************************************************

void htc ( int rank, int np )

//****************************************************************************
{
  int i, fsdo, j, nx, nt, fn, tag, npp;
  double ro, cp, k, t0, tn, x0, xn, h,  dt, fac1, fac2, safe, dtmax;
  double dtlast, time, pos, temp1, temp2, source, sfact, fout, tr, Tsize;
  ofstream t_file;
  double *T, *my_chp, *buddy_chp;
  double *F;
  double *S;
  double *rhs;
  MPI_Status status;

    //initialize values
    x0=0.0;
    xn=100.0;
    h=0.1;
    t0=0.0;
    tn=3.0;
    temp1=0.0;
    temp2=0.0;
    ro=1.0;
    cp=1.0;
    k=1.0;
    safe=0.7;
    sfact=1.0;
    fout=1.0;
    dtmax=1.0;

    //determine number of points
    nx = nearbyint((xn-x0)/h+1);
    //number of points per process
    //must devide exactly 
    npp = nx/np;
    Tsize = npp+4;

    //declare arrays
    T = new double [Tsize];
    my_chp = new double [Tsize];
    buddy_chp = new double [Tsize];
    F = new double [npp+1];
    S = new double [npp];
    rhs = new double [npp];

    //initialize arrays
    for (int i=0; i<Tsize; i++)
        {
            T[i]=0.0;
        }
    for (int i=0; i<npp+1; i++)
        {
            F[i]=0.0;
        }

    //Determining the time step with safety factor
    dt = (h*h/2)*safe;

    //Number of time steps
    nt = floor((tn-t0)/dt);

    //    Applying BC
    if (rank == 0)
      {
	T[0] = temp1;
	T[1] = temp1;
	T[2] = temp1;
      }
    if (rank == np-1)
      {
	T[npp+1] = temp2;
	T[npp+2] = temp2;
        T[npp+3] = temp2;
      }

    fac1 = 15.0/12.0;
    fac2 = 1.0/12.0;

    //
    // Loop over time
    //
    time = 0.0;
    setjmp(restart);
    while (iteration < nt)
    //for ( int j=1; j<nt; j=j+1)
      //for ( int j=1; j<3; j=j+1)
        {
            time = time + dt;
	    //checkpointing every 5 iterations
	    if (0==iteration%5)
	      {
		app_buddy_ckpt(world);
	      }

	    if (0==iteration%51 || rank==np/2)
	      {
		printf("Rank %04d: committing suicide at iteration %d\n", rank, iteration);
		raise(SIGKILL);
              }
		MPI_Barrier(MPI_COMM_WORLD);
	    //
	    // Send T[1],T[2] to rank-1
	    //	    
	    if (rank > 0 && rank < np)
	      {
		tag = 1;
		MPI_Send (&T[2], 2, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
	      }
	    
	    //
	    // Receive T[npp+1],T[npp+2] from rank+1
	    //
	    if (rank < np-1)
	      {
		tag = 1;
		MPI_Recv (&T[npp+2], 2, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
	      }
	    
	    //
	    // Send T[npp-1],T[npp] to rank+1
	    //
	    if (rank < np-1)
	      {
		tag = 1;
		MPI_Send (&T[npp], 2, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
	      }
	    
	    //
	    // Receive T[-2],T[-1] from rank-1
	    //
	    if (rank > 0 && rank < np)
	      {
		tag = 1;
		MPI_Recv (&T[0], 2, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
	      }	    
            
            //Calculating fluxes
	    for (int i=0; i<npp+1; i=i+1)
                {
                    F[i]=fac1*(T[i+2]-T[i+1])-fac2*(T[i+3]-T[i]);
                }

	    //Calculate source
	    double rn=rank*npp;
	    //cout << x0+(2+rn)*h<<"\n";
            for (int i=0; i<npp; i=i+1)
                {
		  pos = x0+(i+rn)*h;
		  //S[i] = 0.01;
		    //                    S[i] = sfact*sin(-(pos-time)*(pos-time));
                    S[i] = sfact*exp(-(pos-time)*(pos-time));
                }
	    //cout << S[0] << "\t" << S[1] << "\t" << S[2] << "\n";
	    //Calculate RHS
            for (int i=0; i<npp; i=i+1)
                {
                    rhs[i] = (k*(F[i+1]-F[i])/(h*h)+S[i])/(ro*cp);
                }
	    //cout << rhs[0] << "\t" << F[1] << "\t" << F[0] << "\n";
	    
	    //Calculate T at t=t+dt
            for (int i=2; i<npp+2; i=i+1)
                {
                    T[i] = T[i]+dt*rhs[i-2];
                }
	    
	    //Reset Boundary values
	    if (rank == 0)
	      {
	    		T[0] = temp1;
	        T[1] = temp1;
	        T[2] = temp1;
	      }
	    if (rank == np-1)
	      {
	    		T[npp+1] = temp2;
	    		T[npp+2] = temp2;
	    	T[npp+3] = temp2;
	      }

	    //cout << j <<"\n";

	    //Output results
	    MPI_Barrier(MPI_COMM_WORLD);
	    char fname[10];
	    sprintf(fname,"%d.%s",100000*(rank+1)+j,"txt");
	    t_file.open (fname);
	    for (int i=2; i<npp+2; i++)
	      {
		t_file << (i-2+rn)*h <<"\t"<< T[i] <<"\t"<< S[i-2] <<"\n";
	      }
	    t_file.close();

	}

  return;
}

//****************************************************************************

// Functions for the FT functionality

//****************************************************************************
//
/* Buddy checkpointing */
//
static int app_buddy_ckpt(MPI_Comm comm) {
    if(0 == rank ) fprintf(stderr, "Rank %04d: checkpointing to %04d after iteration %d\n", rank, rbuddy(rank), iteration);
    /* Store my checkpoint on my "right" neighbor */
    MPI_Sendrecv(T, Tsize, MPI_DOUBLE, rbuddy(rank), ckpt_tag,
                 buddy_chp, Tsize, MPI_DOUBLE, lbuddy(rank), ckpt_tag,
                 comm, MPI_STATUS_IGNORE);
    /* Commit the local changes to the checkpoints only if successful. */
    if(app_needs_repair()) {
        fprintf(stderr, "Rank %04d: checkpoint commit was not succesful, rollback instead\n", rank);
        longjmp(restart, 0);
    }
    ckpt_iteration = iteration;
    /* Memcopy my own memory in my local checkpoint (with datatypes) */
    MPI_Sendrecv(T, Tsize, MPI_DOUBLE, 0, ckpt_tag,
                 my_chp, Tsize, MPI_DOUBLE, 0, ckpt_tag,
                 MPI_COMM_SELF, MPI_STATUS_IGNORE);
    return MPI_SUCCESS;
}
//
/* mockup checkpoint restart: we reset iteration, and we prevent further
 * error injection */
//
static int app_reload_ckpt(MPI_Comm comm) {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);

    /* send my ckpt_iteration to my buddy to decide if we need to exchange a
     * checkpoint
     *   If ckpt_iteration is -1, I am restarting.
     *   It iteration is -1, then my buddy is restarting.
     *   Note: if an error occurs now, it will be absorbed by the error handler
     *   and the restart will be repeated.
     */
    MPI_Sendrecv(&ckpt_iteration, 1, MPI_INT, rbuddy(rank), ckpt_tag,
                 &iteration, 1, MPI_INT, lbuddy(rank), ckpt_tag,
                 comm, MPI_STATUS_IGNORE);

    if( -1 == iteration && -1 == ckpt_iteration ) {
        fprintf(stderr, "Buddy checkpointing cannot restart from this failures because both me and my buddy have lost our checkpoints...\n");
        MPI_Abort(comm, -1);
    }

    if( -1 == iteration ) {
        fprintf(stderr, "Rank %04d: sending checkpoint to %04d at iteration %d\n", rank, lbuddy(rank), ckpt_iteration);
        /* My buddy was dead, send the checkpoint */
        MPI_Send(buddy_chp, Tsize, MPI_DOUBLE, lbuddy(rank), ckpt_tag, comm);
    }
    if( -1 == ckpt_iteration ) {
        /* I replace a dead, get the ckeckpoint */
        fprintf(stderr, "Rank %04d: restarting from %04d at iteration %d\n", rank, rbuddy(rank), iteration);
        MPI_Recv(T, Tsize, MPI_DOUBLE, rbuddy(rank), ckpt_tag, comm, MPI_STATUS_IGNORE);
        /* iteration has already been set by the sendrecv above */
    }
    else {
        /* I am a survivor,
         * Memcopy my own checkpoint back in my memory */
        MPI_Sendrecv(my_chp, Tsize, MPI_DOUBLE, 0, ckpt_tag,
                     T, Tsize, MPI_DOUBLE, 0, ckpt_tag,
                     MPI_COMM_SELF, MPI_STATUS_IGNORE);
        /* Reset iteration */
        iteration = ckpt_iteration;
    }
    return 0;
}

/* repair comm world, reload checkpoints, etc...
 *  Return: true: the app needs to redo some iterations
 *          false: no failure was fixed, we do not need to redo any work.
 */
static int app_needs_repair(void) {
    MPI_Comm tmp;
    MPIX_Comm_replace(world, &tmp);
    if( tmp == world ) return false;
    if( MPI_COMM_NULL != world) MPI_Comm_free(&world);
    world = tmp;
    app_reload_ckpt(world);
    /* Report that world has changed and we need to re-execute */
    return true;
}

/* Do all the magic in the error handler */
static void errhandler_respawn(MPI_Comm* pcomm, int* errcode, ...) {
    int eclass;
    MPI_Error_class(*errcode, &eclass);
    MPI_Error_string(*errcode, estr, &strl);

    if( MPIX_ERR_PROC_FAILED != eclass &&
        MPIX_ERR_REVOKED != eclass ) {
        fprintf(stderr, "%04d: errhandler invoked with unknown error %s\n", rank, estr);
        raise(SIGSEGV);
        MPI_Abort(MPI_COMM_WORLD, *errcode);
    }
    
    fprintf(stderr, "%04d: errhandler invoked with error %s\n", rank, estr);
    
    MPIX_Comm_revoke(*pcomm);
    if(app_needs_repair()) longjmp(restart, 0);
}


static int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm) {
    MPI_Comm icomm, /* the intercomm between the spawnees and the old (shrinked) world */
             scomm, /* the local comm for each sides of icomm */
             mcomm; /* the intracomm, merged from icomm */
    MPI_Group cgrp, sgrp, dgrp;
    int rc, flag, rflag, i, nc, ns, nd, crank, srank, drank;

redo:
    if( comm == MPI_COMM_NULL ) { /* am I a new process? */
        /* I am a new spawnee, waiting for my new rank assignment
         * it will be sent by rank 0 in the old world */
        MPI_Comm_get_parent(&icomm);
        scomm = MPI_COMM_WORLD;
        MPI_Recv(&crank, 1, MPI_INT, 0, 1, icomm, MPI_STATUS_IGNORE);
        if( verbose ) {
            MPI_Comm_rank(scomm, &srank);
            fprintf(stderr, "Spawnee %d: crank=%d\n", srank, crank);
        }
    }
    else {
        /* I am a survivor: Spawn the appropriate number
         * of replacement processes (we check that this operation worked
         * before we procees further) */
        /* First: remove dead processes */
        MPIX_Comm_shrink(comm, &scomm);
        MPI_Comm_size(scomm, &ns);
        MPI_Comm_size(comm, &nc);
        nd = nc-ns; /* number of deads */
        if( 0 == nd ) {
            /* Nobody was dead to start with. We are done here */
            MPI_Comm_free(&scomm);
            *newcomm = comm;
            return MPI_SUCCESS;
        }
        /* We handle failures during this function ourselves... */
        MPI_Comm_set_errhandler( scomm, MPI_ERRORS_RETURN );

        rc = MPI_Comm_spawn(gargv[0], &gargv[1], nd, MPI_INFO_NULL,
                            0, scomm, &icomm, MPI_ERRCODES_IGNORE);
        flag = (MPI_SUCCESS == rc);
        MPIX_Comm_agree(scomm, &flag);
        if( !flag ) {
            if( MPI_SUCCESS == rc ) {
                MPIX_Comm_revoke(icomm);
                MPI_Comm_free(&icomm);
            }
            MPI_Comm_free(&scomm);
            if( verbose ) fprintf(stderr, "%04d: comm_spawn failed, redo\n", rank);
            goto redo;
        }

        /* remembering the former rank: we will reassign the same
         * ranks in the new world. */
        MPI_Comm_rank(comm, &crank);
        MPI_Comm_rank(scomm, &srank);
        /* the rank 0 in the scomm comm is going to determine the
         * ranks at which the spares need to be inserted. */
        if(0 == srank) {
            /* getting the group of dead processes:
             *   those in comm, but not in scomm are the deads */
            MPI_Comm_group(comm, &cgrp);
            MPI_Comm_group(scomm, &sgrp);
            MPI_Group_difference(cgrp, sgrp, &dgrp);
            /* Computing the rank assignment for the newly inserted spares */
            for(i=0; i<nd; i++) {
                MPI_Group_translate_ranks(dgrp, 1, &i, cgrp, &drank);
                /* sending their new assignment to all new procs */
                MPI_Send(&drank, 1, MPI_INT, i, 1, icomm);
            }
            MPI_Group_free(&cgrp); MPI_Group_free(&sgrp); MPI_Group_free(&dgrp);
        }
    }

    /* Merge the intercomm, to reconstruct an intracomm (we check
     * that this operation worked before we proceed further) */
    rc = MPI_Intercomm_merge(icomm, 1, &mcomm);
    rflag = flag = (MPI_SUCCESS==rc);
    MPIX_Comm_agree(scomm, &flag);
    if( MPI_COMM_WORLD != scomm ) MPI_Comm_free(&scomm);
    MPIX_Comm_agree(icomm, &rflag);
    MPI_Comm_free(&icomm);
    if( !(flag && rflag) ) {
        if( MPI_SUCCESS == rc ) {
            MPI_Comm_free(&mcomm);
        }
        if( verbose ) fprintf(stderr, "%04d: Intercomm_merge failed, redo\n", rank);
        goto redo;
    }

    /* Now, reorder mcomm according to original rank ordering in comm
     * Split does the magic: removing spare processes and reordering ranks
     * so that all surviving processes remain at their former place */
    rc = MPI_Comm_split(mcomm, 1, crank, newcomm);

    /* Split or some of the communications above may have failed if
     * new failures have disrupted the process: we need to
     * make sure we succeeded at all ranks, or retry until it works. */
    flag = (MPI_SUCCESS==rc);
    MPIX_Comm_agree(mcomm, &flag);
    MPI_Comm_free(&mcomm);
    if( !flag ) {
        if( MPI_SUCCESS == rc ) {
            MPI_Comm_free( newcomm );
        }
        if( verbose ) fprintf(stderr, "%04d: comm_split failed, redo\n", rank);
        goto redo;
    }

    /* restore the error handler */
    if( MPI_COMM_NULL != comm ) {
        MPI_Errhandler errh;
        MPI_Comm_get_errhandler( comm, &errh );
        MPI_Comm_set_errhandler( *newcomm, errh );
    }

    return MPI_SUCCESS;
}
