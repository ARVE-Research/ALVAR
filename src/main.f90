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
! $ mpirun -np 18 ./src/alvar /home/terraces/datasets/dgvm_input/climate/transient1871-2010_list-formatted.nc -110.25/40.25 1990/10 test.nc
!
! NOTE: user-specified lon/lat must be **EXACT of the input file values** (routine to find closest gridcell is yet to be written) - Leo Lai 28 Apr 2021


! Main program to distribute work among cores

use parametersmod, only : i1,i2,i4,sp,dp
use errormod,      only : ncstat,netcdf_err
use coordsmod,     only : index,parsecoords
use outputmod,     only : infompi,getoutfile
use modelmod,      only : model
use gwgenmod,      only : gwgen
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

type(infompi) :: info

integer(i4), allocatable, dimension(:) :: srt
integer(i4), allocatable, dimension(:) :: cnt

integer(i1), allocatable, dimension(:) :: ob

integer(i4), dimension(2) :: job

!--------------------

call MPI_INIT(ierr)

call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

!--------------------

allocate(srt(numtasks))
allocate(cnt(numtasks))

srt = 0
cnt = 0

info%nproc = numtasks

infosize = sizeof(info)

sendsize = infosize + (8*numtasks)

allocate(ob(sendsize))

ob = 0

!--------------------

if (rank == 0) then

  call startmpi(info, srt, cnt)

  info%validcell = sum(cnt)

  call getoutfile(info)

  call infotobyte(info, srt, cnt, ob)

end if

!---
! Broadcast all info to all processes

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

call MPI_BCAST(ob, (sendsize/4), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

!---

if (rank == 0) write(0,*) 'Broadcast complete'

if (rank /= 0) then

  call bytetoinfo(info, srt, cnt, ob)

end if

!---

job = [srt(rank+1), cnt(rank+1)]

write(0,*) 'Rank:', rank, 'recieved cell srt and cnt: ', job

!--------------------

call model(info, job, rank)

!--------------------

call MPI_FINALIZE(ierr)



!-------------------------------------------------------

contains

subroutine startmpi(info,srt,cnt)

use statevarsmod, only : startyr,calcyrs

type(infompi), target              , intent(inout) :: info
integer(i4),           dimension(:), intent(inout) :: srt
integer(i4),           dimension(:), intent(inout) :: cnt

integer :: i
integer :: n

integer :: ifid
integer :: dimid
integer :: ilen
integer :: tlen

!--- File variables
character(100), pointer :: infile
character(100), pointer :: outfile
character(100), pointer :: timestring
integer(i4)   , pointer :: nproc
integer(i4)   , pointer :: t0
integer(i4)   , pointer :: nt
real(sp)      , pointer :: plon
real(sp)      , pointer :: plat

character(100) :: gridlonlat

integer :: len
integer :: remainder

type(index) :: timevals

integer(i4) :: baseyr
integer(i4) :: t1


infile     => info%infile
outfile    => info%outfile
timestring => info%timestring
nproc      => info%nproc
t0         => info%t0
nt         => info%nt
plon       => info%plon
plat       => info%plat


call getarg(1,infile)

ncstat = nf90_open(infile,nf90_nowrite,ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'index',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ilen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!--------------------

call getarg(2,gridlonlat)

call parsecoords(gridlonlat,timevals)

plon = real(timevals%minlon)
plat = real(timevals%minlat)

if (plon > 180. .OR. plon < -180) then
  write(0,*) ' '
  write(0,*) 'User-specified longitude exceeded dimension'
  write(0,*) ' '
  stop
end if

if (plat > 90. .OR. plat < -90) then
  write(0,*) ' '
  write(0,*) 'User-specified latitude exceeded dimension'
  write(0,*) ' '
  stop
end if

!--------------------

baseyr = 1871

call getarg(3,timestring)

call parsecoords(timestring,timevals)

startyr = nint(timevals%minlon)
calcyrs = nint(timevals%minlat)

t0 = 1 + 12 * (startyr - baseyr)
t1 = t0 + 12 * calcyrs - 1

nt = t1 - t0 + 1

!--------------------

len = floor(real(ilen) / nproc)

remainder = ilen - (len * nproc)

n = 1
do i = 1, nproc

  srt(n) = (i-1) * len + 1
  cnt(n) = len

  if (i == nproc) then

    cnt(n) = len + remainder

  end if

  n = n + 1

end do

!--------------------

call getarg(4, outfile)


end subroutine startmpi

!-------------------------------------------------------

subroutine infotobyte(info,srt,cnt,ob)

type(infompi)              , intent(in)    :: info
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

type(infompi)              , intent(inout)    :: info
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
