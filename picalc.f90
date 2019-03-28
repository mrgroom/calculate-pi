! Program to calculate pi using a parallel Monte Carlo method
! to compile (using intel fortran): mpifort -O3 -fp-model precise -o CalulatePi picalc.f90
program CalculatePi
use mpi_f08
implicit none
integer :: N,i,c
real(8) :: xc,yc,r,pi,rPi,rDummy,rSum
real(8),dimension(2) :: rnum
integer :: iStart,iEnd,iOffset,iSize
! MPI variables
integer :: iRank, iNumTasks, iErr, iTag, status(MPI_STATUS_SIZE) 

! Initialise MPI
call MPI_INIT(iErr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, iNumTasks, iErr)
call MPI_COMM_RANK(MPI_COMM_WORLD, iRank, iErr)

! Problem parameters
N = 100000000 ! NUmber of points
r = 1.d0 ! Radius of circle

! Make sure problem is distributed equally
if(mod(N,iNumTasks).ne.0.and.iRank.eq.0) then
   print*,'Number of points N =',N,' is not divisible by number of mpi processes n =',iNumTasks     
   call MPI_ABORT(MPI_COMM_WORLD,1,iErr)
   stop
end if
call MPI_BARRIER(MPI_COMM_WORLD)

! Decompose loop
iSize = N/iNumTasks
iStart = 1 + iRank*iSize
iEnd = (iRank+1)*iSize
!print*,'Rank =',iRank,'Start =',iStart,'End =',iEnd

! Calculate pi
call random_seed()
c = 0
do i=iStart,iEnd
    call random_number(rnum)
    xc = -r+rnum(1)*2.d0*r
    yc = -r+rnum(2)*2.d0*r
    if(xc**2+yc**2.le.r**2) c = c+1
end do

! Collect results
rDummy = 0.d0
call MPI_ALLREDUCE(dble(c),rDummy,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iErr)
rSum = rDummy

! Final answer
pi = 4.d0*rSum/real(N,kind=8)
rPi = 4.d0*atan(1.d0) ! System value of Pi
if(iRank.eq.0) then
   print*,'Calculated value of Pi is: ',pi
   print*,'Error compared to system value of pi: ',abs(pi-rPi)
end if
call MPI_BARRIER(MPI_COMM_WORLD)

! Finalise MPI
call MPI_FINALIZE(iErr)

end program CalculatePi
