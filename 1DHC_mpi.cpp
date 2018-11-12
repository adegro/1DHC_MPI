# include <cmath>
# include <cstdlib>
# include <ctime>
# include <fstream>
# include <iostream>
# include <mpi.h>
# include <mpi-ext.h>
# include <signal.h>

using namespace std;

int main ( int argc, char *argv[] );
void htc ( int rank, int np );
static void verbose_errhandler(MPI_Comm* pcomm, int* perr, ...);
static int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm);
static int app_buddy_ckpt(MPI_Comm comm);
static int app_reload_ckpt(MPI_Comm comm);
static int app_needs_repair(void);

//****************************************************************************

int main ( int argc, char *argv[] )

//****************************************************************************
{
  int rank;
  int size, np;
  MPI_Errhandler errh;

//
//  Initialize MPI.
//

  MPI_Init ( &argc, &argv );
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );
  MPI_Comm_create_errhandler(verbose_errhandler,
                               &errh);
  MPI_Comm_set_errhandler(MPI_COMM_WORLD,
                            errh);

  np = size;
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
  double dtlast, time, pos, temp1, temp2, source, sfact, fout, tr;
  ofstream t_file;
  double *T;
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

    //declare arrays
    T = new double [npp+4];
    F = new double [npp+1];
    S = new double [npp];
    rhs = new double [npp];

    //initialize arrays
    for (int i=0; i<npp+4; i++)
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
    for ( int j=1; j<nt; j=j+1)
      //for ( int j=1; j<3; j=j+1)
        {
            time = time + dt;

	    if (j == 50)
	      {

		MPI_Barrier(MPI_COMM_WORLD);
		if( rank == (np-1)||rank == (np/2) )
		  {
		    printf("Rank %d / %d: bye bye!\n", rank, np);
		    raise(SIGKILL);
		  }

		MPI_Barrier(MPI_COMM_WORLD);
		printf("Rank %d / %d: Stayin' alive!\n", rank, np);
              }
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

static void verbose_errhandler(MPI_Comm* pcomm, int* perr, ...) {

//****************************************************************************
    MPI_Comm comm = *pcomm;
    int err = *perr;
    char errstr[MPI_MAX_ERROR_STRING];
    int i, rank, size, nf, len, eclass;
    MPI_Group group_c, group_f;
    int *ranks_gc, *ranks_gf;

    MPI_Error_class(err, &eclass);
    if( MPIX_ERR_PROC_FAILED != eclass ) {
        MPI_Abort(comm, err);
    }

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* We use a combination of 'ack/get_acked' to obtain the list of 
     * failed processes (as seen by the local rank). 
     */
    MPIX_Comm_failure_ack(comm);
    MPIX_Comm_failure_get_acked(comm, &group_f);
    MPI_Group_size(group_f, &nf);
    MPI_Error_string(err, errstr, &len);
    printf("Rank %d / %d: Notified of error %s. %d found dead: { ",
           rank, size, errstr, nf);

    /* We use 'translate_ranks' to obtain the ranks of failed procs 
     * in the input communicator 'comm'.
     */
    ranks_gf = (int*)malloc(nf * sizeof(int));
    ranks_gc = (int*)malloc(nf * sizeof(int));
    MPI_Comm_group(comm, &group_c);
    for(i = 0; i < nf; i++)
        ranks_gf[i] = i;
    MPI_Group_translate_ranks(group_f, nf, ranks_gf,
                              group_c, ranks_gc);
    for(i = 0; i < nf; i++)
        printf("%d ", ranks_gc[i]);
    printf("}\n");
    free(ranks_gf); free(ranks_gc);
}
//
/* Buddy checkpointing */
//
static int app_buddy_ckpt(MPI_Comm comm) {
    if(0 == rank || verbose) fprintf(stderr, "Rank %04d: checkpointing to %04d after iteration %d\n", rank, rbuddy(rank), iteration);
    /* Store my checkpoint on my "right" neighbor */
    MPI_Sendrecv(mydata_array, count, MPI_DOUBLE, rbuddy(rank), ckpt_tag,
                 buddy_ckpt,   count, MPI_DOUBLE, lbuddy(rank), ckpt_tag,
                 comm, MPI_STATUS_IGNORE);
    /* Commit the local changes to the checkpoints only if successful. */
    if(app_needs_repair()) {
        fprintf(stderr, "Rank %04d: checkpoint commit was not succesful, rollback instead\n", rank);
        longjmp(restart, 0);
    }
    ckpt_iteration = iteration;
    /* Memcopy my own memory in my local checkpoint (with datatypes) */
    MPI_Sendrecv(mydata_array, count, MPI_DOUBLE, 0, ckpt_tag,
                 my_ckpt, count, MPI_DOUBLE, 0, ckpt_tag,
                 MPI_COMM_SELF, MPI_STATUS_IGNORE);
    return MPI_SUCCESS;
}

