program main

! Before compiling
! $ module load netCDF-Fortran/4.5.3-gmpich-2021.01 Autotools pkgconfig

! To compile:
! $ ./autogen.sh
! $ ./configure
! $ make

! To run: program was written to run ALL validcells globally, but number of years can be specified
! $ module purge              ! because Autotools switched GCC back to 9.3.0. v10.2 is required to run the code
! $ module load netCDF-Fortran/4.5.3-gmpich-2021.01
!
! $ mpirun -np 48 ./src/alvar /home/terraces/datasets/dgvm_input/climate/transient1871-2010_list-formatted.nc -110.25/40.25 1990/5 test.nc
!
! NOTE: user-specified lon/lat must be **EXACT of the input file values** (routine to find closest gridcell is yet to be written) - Leo Lai 28 Apr 2021


! Main program to distribute work among cores

use parametersmod, only : i1,i2,i4,sp,dp
use errormod,      only : ncstat,netcdf_err
use coordsmod,     only : index,parsecoords
use mpivarsmod,    only : mpivars
use initmod,       only : initjob
use outputmod,     only : getoutfile,getoutfile_onelayer
use drivermod,     only : driver
use netcdf
use mpi

implicit none

!--------------------
! Parallel program to read in list-formatted ncfile an calculate the average of each validcell across specified timeframe
integer :: rank
integer :: numtasks
integer :: ierr

integer :: infosize
integer :: sendsize

type(mpivars) :: info

integer(i4), allocatable, dimension(:) :: srt
integer(i4), allocatable, dimension(:) :: cnt

integer(i1), allocatable, dimension(:) :: ob

integer(i4), dimension(2) :: job

!--------------------

call MPI_INIT(ierr)

call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

!--------------------

allocate(srt(numtasks))
allocate(cnt(numtasks))

srt = 0
cnt = 0

info%nproc = numtasks

infosize = sizeof(info)

sendsize = infosize + (8*numtasks)

!--------------------

allocate(ob(sendsize))

ob = 0

!--------------------

if (rank == 0) then

  call initjob(info,srt,cnt)

  ! call getoutfile(info)

  call getoutfile_onelayer(info)

  call infotobyte(info,srt,cnt,ob)

end if

!--------------------
! Broadcast all info to all processes

call mpi_barrier(MPI_COMM_WORLD,ierr)

call mpi_bcast(ob,(sendsize/4),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!--------------------

if (rank == 0) write(0,*) 'Broadcast complete'

if (rank /= 0) then

  call bytetoinfo(info,srt,cnt,ob)

end if

!--------------------

job = [srt(rank+1),cnt(rank+1)]

write(0,*) 'Rank:', rank, 'recieved cell srt and cnt: ', job

!--------------------

call driver(info,job,rank)

!--------------------

call mpi_finalize(ierr)

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

contains

!-------------------------------------------------------

subroutine infotobyte(info,srt,cnt,ob)

type(mpivars)              , intent(in)    :: info
integer(i4)  , dimension(:), intent(in)    :: srt
integer(i4)  , dimension(:), intent(in)    :: cnt
integer(i1)  , dimension(:), intent(inout) :: ob

integer(i4) :: len
integer(i4) :: ilen
integer(i4) :: alen
integer(i4) :: n
integer(i4) :: m

ilen = sizeof(info)

alen = sizeof(srt)

len = sizeof(ob)

n = ilen

ob(1:n) = transfer(info, ob(1:n))

m = n + alen

ob((n+1):m) = transfer(srt, ob((n+1):m))

n = m + alen

ob((m+1):n) = transfer(cnt, ob((m+1):n))

end subroutine infotobyte

!-------------------------------------------------------

subroutine bytetoinfo(info,srt,cnt,ob)

type(mpivars)              , intent(inout)    :: info
integer(i4)  , dimension(:), intent(inout)    :: srt
integer(i4)  , dimension(:), intent(inout)    :: cnt
integer(i1)  , dimension(:), intent(in)       :: ob

integer :: len
integer :: ilen
integer :: alen
integer :: n
integer :: m

ilen = sizeof(info)

alen = sizeof(srt)

len = sizeof(ob)

n = ilen

info = transfer(ob(1:n), info)

m = n + alen

srt = transfer(ob((n+1):m), srt)

n = m + alen

cnt = transfer(ob((m+1):n), cnt)

end subroutine bytetoinfo

!-------------------------------------------------------

end program main
