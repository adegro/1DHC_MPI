!--------------------------------------------------
program fdhc

  include 'mpif.h'
  include 'mpif-ext.h'
!  implicit none
  
  integer :: rank, size, np, ierr, errstr(MPI_MAX_ERROR_STRING)
  integer :: i, fsdo, j, nx, nt, fn, npp
  integer :: stat(MPI_STATUS_SIZE), tag
  real*8 :: ro, cp, k, t0, tn, x0, xn, h,  dt, fac1, fac2, safe
  real*8 :: dtlast, time, pos, temp1, temp2, source, sfact, rn
  real*8,dimension(:),allocatable :: T, F, S, rhs
  character(LEN=3) :: tinteg
  integer :: errh

  !
  !Initialize MPI
  !
  call MPI_INIT (ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size,  ierr)
  call MPI_Comm_create_errhandler(verbose_errhandler, &errh, ierr);
  call MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN, ierr)

  np = size

  !     read input parameters
  call inpar(x0,xn,h,t0,tn,temp1,temp2, &
       ro,cp,k,safe,sfact, &
       fsdo,tinteg)

  !Determine number of points
  nx = nint((xn-x0)/h+1)
  !Number of points per process
  !must devide exactly
  npp = nx/np

  !Allocating memory
  allocate (T(npp+4))
  allocate (F(npp+1))
  allocate (S(npp))
  allocate (rhs(npp))

  !Initializing arrays
  T = 0
  F = 0

  !Determining the time step with safety factor
  dt = (h*h/2)*safe

  !Number of time steps
  nt = floor((tn-t0)/dt)

  !Applying BC
  if (rank==0) then
     T(1) = temp1
     T(2) = temp1
     T(3) = temp1
  else if (rank==np-1) then
     T(npp+2) = temp2
     T(npp+3) = temp2
     T(npp+4) = temp2
  end if

  fac1 = 15.0d+00/12.0d+00
  fac2 = 1.0d+00/12.0d+00

  !----------------------------------------
  ! Main loop over time
  !----------------------------------------

  time = 0

  do j = 1,nt
  !do j=1,3

     time = time + dt

     !!
     !!Exchange information between domains if size>1
     !!
     !
     !Send T[1],T[2] to rank-1
     !
     if (size == 1) then
        goto 1000
     end if
     
     if (rank > 0 .and. rank < np) then
        tag = 1
        call MPI_Send (T(3), 2, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, ierr)
     end if
     
     !
     !Receive T[npp+1],T[npp+2] from rank+1
     !
     if (rank < np-1) then
        tag = 1
        call MPI_Recv (T(npp+3), 2, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, stat, ierr)
     end if
     
     !Send T[npp-1],T[npp] to rank+1
     !
     if (rank < np-1) then
        tag = 1
        call MPI_Send (T(npp+1), 2, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, ierr)
     end if

     !
     !Receive T[-2],T[-1] from rank-1
     !
     if (rank > 0 .and. rank < np) then
        tag = 1
        call MPI_Recv (T(1), 2, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, stat, ierr)
     end if

     1000 continue
     
     !
     !Calculating fluxes
     !
     do i = 2,npp+2
     !do i = 2,8
        call flux_4tho(F(i-1),T(i+1),T(i),T(i-1),T(i+2),fac1,fac2)
      !  write(6,*) F(i-1),T(i+1),T(i),T(i-1),T(i+2),fac1,fac2
     end do

     !Calculating source
     rn=rank*npp
     do i = 2,npp+1
        pos = x0+((i-2.0d+00)+rn)*h
        S(i-1) = source(pos,time,sfact)
        !S(i-1) = 0.1
     end do

     !Calculating RHS
     do i = 3,npp+2
        rhs(i-2) = (k*(F(i-1)-F(i-2))/(h*h)+S(i-2))/(ro*cp)
     end do

     !Euler time integration
     do i = 3,npp+2
        T(i)=T(i)+dt*rhs(i-2)
     end do

     !Write out results
     fn=100000*rank+j
     call output(time,T,npp,fn,fsdo,h,rn)

  end do
end program fdhc
!
!------------------------------------------------------------------
!******************************************************************
!------------------------------------------------------------------
!Functions and subroutines
!------------------------------------------------------------------
!******************************************************************
!------------------------------------------------------------------
  !------------------------------------------------------------------
  real*8 function source(x,t,sfact)
    real*8 :: x,t,sfact
    source = sfact*exp(-(x-t)*(x-t))
    return
  end function source
  !------------------------------------------------------------------
  subroutine flux_4tho(f,tp1,t1,tm1,tp2,fac1,fac2)
    real*8 :: f, tp1, t1, tm1, tp2, fac1, fac2
    f=fac1*(tp1-t1)-fac2*(tp2-tm1)
    return
  end subroutine flux_4tho
  !------------------------------------------------------------------
  subroutine output(time,T,size,fn,fsdo,h,rn)
    integer :: fn,fsdo,size
    real*8 :: time,h,rn
    real*8, dimension(size+2) :: T
    character(LEN=10) :: fname

    !file name
    write(fname,200) fn,'.out'

    !file number
    fn = 10000+fn

    !open file
    open(unit = fn, file = fname)

    do i=2,size+1
       write(fn,100) (i-2+rn)*h,T(i)
    end do

    close(fn)

100 format(2f9.4)
200 format(i0.6,a)
  end subroutine output
  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine inpar(x0,xn,h,t0,tn,temp1,temp2, &
       ro,cp,k,safe,sfact, &
       fsdo,tinteg)

    implicit none

    integer :: nlines, iline, fsdo
    real*8 :: x0,xn,h,t0,tn,temp1,temp2
    real*8 :: ro,cp,k,safe,sfact
    real*8 :: text(20)
    character(len=3) :: tinteg

    open(unit=10,file='inp.par',err=101,status='old')
    rewind 10

    !
    !     -----this sub reads in the input parameters for the solver
    !
    read(10,*)nlines

    do iline=1,nlines
       read(10,1) text
    end do

    read(10,1) text
    read(10,*) x0,xn,h,t0,tn,temp1,temp2

    read(10,1) text
    read(10,*) ro,cp,k,safe,sfact

    read(10,1) text
    read(10,*) fsdo, tinteg

1   format(20a4)

    !
    !     -----as we found all the files, and were able to open them:
    !
    goto 201

101 continue
    write(6,*)' could not find input file'
    !      exit solver
    stop

201 continue
    return
  end subroutine inpar
  !------------------------------------------------------------------
