module initmod

implicit none

contains

!---------------------------------------------------------------------

subroutine initjob(info,srt,cnt)

use parametersmod, only : i4
use errormod,      only : ncstat,netcdf_err
use coordsmod,     only : index,parsecoords
use mpivarsmod,    only : mpivars
use statevarsmod,  only : calcyrs
use netcdf

implicit none

type(mpivars), target,               intent(inout) :: info
integer(i4),           dimension(:), intent(inout) :: srt
integer(i4),           dimension(:), intent(inout) :: cnt

character(200), pointer :: outfile
character(200), pointer :: gridlistfile
character(200), pointer :: cfile_spinup
character(200), pointer :: cfile_transient
character(200), pointer :: soilfile
character(200), pointer :: topofile
logical,        pointer :: dospinup
logical,        pointer :: dotransient
integer(i4),    pointer :: spin_baseyr
integer(i4),    pointer :: spin_startyr
integer(i4),    pointer :: spinupyears
integer(i4),    pointer :: tran_baseyr
integer(i4),    pointer :: tran_startyr
integer(i4),    pointer :: transientyears
integer(i4),    pointer :: ilen
integer(i4),    pointer :: spin_tlen
integer(i4),    pointer :: tran_tlen
character(200), pointer :: timestring
integer(i4)   , pointer :: outmode
integer(i4)   , pointer :: nproc

character(200) :: jobfile
integer :: ifid
integer :: dimid

character(1) :: outmode_str
integer :: len
integer :: remainder

type(index) :: timevals

integer(i4) :: t1

integer :: i
integer :: n



namelist /joboptions/ &
  cfile_spinup,       &
  cfile_transient,    &
  soilfile,           &
  topofile,           &
  dospinup,           &
  dotransient,        &
  spin_baseyr,        &
  spin_startyr,       &
  spinupyears,        &
  tran_baseyr
!   fixedco2,           &
!   ocean_uptake,       &
!   cal_year,           &
!   co2file,            &
!   popdfile,           &
!   poppfile,           &
!   pftparsfile,        &
!   nspinyrsout,        &
!   outputvar,          &
!   lu_turn_yrs,        &
!   maxmem,             &
!   nolanduse,          &
!   startyr_foragers

!--------------------

outfile         => info%outfile
gridlistfile    => info%gridlistfile
cfile_spinup    => info%cfile_spinup
cfile_transient => info%cfile_transient
soilfile        => info%soilfile
topofile        => info%topofile
dospinup        => info%dospinup
dotransient     => info%dotransient
spin_baseyr     => info%spin_baseyr
spin_startyr    => info%spin_startyr
spinupyears     => info%spinupyears
tran_baseyr     => info%tran_baseyr
tran_startyr    => info%tran_startyr
transientyears  => info%transientyears
ilen            => info%ilen
spin_tlen       => info%spin_tlen
tran_tlen       => info%tran_tlen

timestring      => info%timestring
outmode         => info%outmode
nproc           => info%nproc

!--------------------
! Initialize variables with a default value if they are not specified in the namelist
spinupyears    = -9999
transientyears = 1
! nspinyrsout    = -9999
! nolanduse      = .false.
! startyr_foragers = 1000

!--------------------

call getarg(1,jobfile)

open(21,file=jobfile,status='old')

read(21,nml=joboptions)

close(21)

!--------------------

ncstat = nf90_open(cfile_spinup,nf90_nowrite,ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'index',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ilen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!--------------------
! File containing list of grid lon/lat for single cell runs

call getarg(2,gridlistfile)

!--------------------
! Output mode of the model, 0 for full grid array and 1 for single grids

call getarg(3,outmode_str)

read(outmode_str,'(I1)') outmode

!--------------------
! Timestring to indicate start and calcyears for the transient run

call getarg(4,timestring)

call parsecoords(timestring,timevals)

tran_startyr   = nint(timevals%minlon)
transientyears = nint(timevals%minlat)

!--------------------
! Output netCDF file

call getarg(5,outfile)

!--------------------
! Calculate the start and count of gridcell for each MPI process

len = floor(real(ilen) / nproc)

remainder = ilen - (len * nproc)

srt(1)           = 1
cnt              = len
cnt(1:remainder) = len + 1

do i = 2, nproc

  srt(i) = srt(i-1) + cnt(i-1)

end do

end subroutine initjob

!---------------------------------------------------------------------

subroutine initdate(infile,rank,baseyr,startyr,calcyrs,tlen,t0,p0,p1,nt,cntt,nd)

! Initiate time-related variables for model
! start & count of dimensions, startyr, endyr, calcyr and nd (number of days in each month)

use parametersmod, only : i4,ndaymonth
use errormod,      only : ncstat,netcdf_err
use netcdf
use mpi

implicit none

character(200), intent(in)    :: infile
integer(i4),    intent(in)    :: rank
integer(i4),    intent(in)    :: baseyr
integer(i4),    intent(in)    :: startyr
integer(i4),    intent(in)    :: calcyrs
integer(i4),    intent(out)   :: tlen
integer(i4),    intent(out)   :: t0
integer(i4),    intent(out)   :: p0
integer(i4),    intent(out)   :: p1
integer(i4),    intent(out)   :: nt
integer(i4),    intent(out)   :: cntt
integer(i4),    allocatable, dimension(:), intent(inout) :: nd

! Local variables
integer(i4) :: endyr
integer(i4) :: t1
integer(i4) :: yr
integer(i4) :: mon
integer(i4) :: i

integer(i4) :: ifid
integer(i4) :: dimid
integer(i4) :: varid

!----------------------------------------------------

ncstat = nf90_open(infile,nf90_nowrite,ifid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           ! Open netCDF-file (inpput file name, no writing rights, assigned file number)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! Check for errors (after every step)

ncstat = nf90_inq_dimid(ifid,'time',dimid)             ! Get dimension ID for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)   ! Get length of dimension 'time' and assign it to variable tlen
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!----------------------------------------------------
! Get the buffered start and end of the time dimension of each validpixel

t0 = 1 + 12 * (startyr - baseyr)
t1 = t0 + 12 * calcyrs - 1

nt = t1 - t0 + 1

cntt = 12 * (calcyrs + 2)          ! months, includes one-year (12 months) buffer on either end

!-------------------
! calculate file start and count indices

t1 = t0 + 12 * calcyrs - 1         ! index of the last month

if (t0 == 1) then                  ! there is no additional data to be had at the front end so copy the first year twice
  p0 = 13
else
  p0 = 1                           ! there is additional data before first year so grab it
  t0 = t0 - 12
end if

if (t1 == tlen) then               ! there is no additional data to be had at the back end so copy the last year twice
  p1 = cntt - 12
else
  p1 = cntt
  t1 = t1 + 12
end if

nt = t1 - t0 + 1          ! nt changes to actual count of number of elements including buffer years (+/- 1 year)

!-------------------

if (rank == 0) then
  write(0,*) startyr,calcyrs
  write(0,*) cntt,nt
  write(0,*) t0,t1
  write(0,*) p0,p1
end if

!-------------------
! Calculate days per month

endyr = startyr + calcyrs - 1

allocate(nd(cntt))

i = 1
do yr = startyr-1, endyr+1
  do mon = 1, 12

    nd(i) = ndaymonth(yr,mon)

    i = i + 1

  end do
end do

end subroutine initdate

!---------------------------------------------------------------------

subroutine initinvars(gridcount,cntt,invars)

use parametersmod, only : i4
use statevarsmod,  only : in_metvars

implicit none

integer(i4),                                 intent(in)    :: gridcount
integer(i4),                                 intent(in)    :: cntt
type(in_metvars), allocatable, dimension(:), intent(inout) :: invars

integer :: i

!-------------------

allocate(invars(gridcount))

do i = 1, gridcount

  allocate(invars(i)%tmp(cntt))
  allocate(invars(i)%dtr(cntt))
  allocate(invars(i)%pre(cntt))
  allocate(invars(i)%wet(cntt))
  allocate(invars(i)%cld(cntt))
  allocate(invars(i)%wnd(cntt))

end do


end subroutine initinvars

!---------------------------------------------------------------------

subroutine initlonlat(infile,gridstart,gridcount,xlen,ylen,inlon,inlat,indx,lon,lat)

! 1. Save current lon and lat of the gridcells into statevars lon and lat
! 2. Decide if the current CPU process would require to print

use parametersmod, only : i4,sp,dp
use errormod,      only : ncstat,netcdf_err
use netcdf
use mpi

implicit none

character(200),                              intent(in)  :: infile
integer(i4),                                 intent(in)  :: gridstart
integer(i4),                                 intent(in)  :: gridcount
integer(i4),                                 intent(out) :: xlen
integer(i4),                                 intent(out) :: ylen
real(dp),       allocatable, dimension(:),   intent(out) :: inlon
real(dp),       allocatable, dimension(:),   intent(out) :: inlat
integer(i4),    allocatable, dimension(:,:), intent(out) :: indx
real(dp),                    dimension(:),   intent(out) :: lon
real(dp),                    dimension(:),   intent(out) :: lat

! Local variables
integer(i4) :: grid
integer(i4) :: ll
integer(i4) :: lon_loc
integer(i4) :: lat_loc
integer(i4), dimension(2) :: ll_loc

integer(i4) :: ifid
integer(i4) :: dimid
integer(i4) :: varid
integer(i4) :: i

! -----------------------------------------------------
! INPUT: Read dimension IDs and lengths of dimensions

ncstat = nf90_open(infile,nf90_nowrite,ifid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           ! Open netCDF-file (inpput file name, no writing rights, assigned file number)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! Check for errors (after every step)

ncstat = nf90_inq_dimid(ifid,'lon',dimid)              ! get dimension ID from dimension 'lon' in the input file
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lat',dimid)              ! get dimension ID from dimension 'lon' in the input file
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----------------------------------------------------
! Read variable IDs and values

allocate(inlon(xlen))       ! Allocate length to longitude array
allocate(inlat(ylen))       ! Allocate length to latitude array
allocate(indx(xlen,ylen))   ! Allocate length to index array

ncstat = nf90_inq_varid(ifid,"lon",varid)                ! Get variable ID for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,inlon)                  ! Get variable values for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)                ! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,inlat)                  ! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"index",varid)              ! Get variable ID for index
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,indx)                   ! Get variable values for index
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---------------------

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---------------------

do grid = 1, gridcount

  ll = gridstart + (grid - 1)   ! Get the current index value

  ! Get the value of lon and lat from index dimension
  ll_loc = findloc(indx, ll)

  lon_loc = ll_loc(1)
  lat_loc = ll_loc(2)

  lon(grid) = inlon(lon_loc)
  lat(grid) = inlat(lat_loc)

end do

end subroutine initlonlat

!---------------------------------------------------------------------

subroutine calcndyear(year,nd,ndyear)

! Find number of days in current year (365 or 366)

use parametersmod, only : i4,sp

implicit none

integer(i4),               intent(in)  :: year
integer(i4), dimension(:), intent(in)  :: nd
integer(i4),               intent(out) :: ndyear

integer(i4) :: start
integer(i4) :: end

!-------------------

start = (12 * year) + 1     ! Start month (i.e. Jan) begin after 12 months of the buffer year
end   = start + 11          ! End month (i.e. Dec)

!-------------------
! Calculate ndyear from sum of nd
ndyear  = sum(nd(start:end))

end subroutine calcndyear

!---------------------------------------------------------------------

subroutine copymonvars(year,grid,tmp,dtr,pre,wet,cld,wnd,ynd)

use parametersmod, only : i4,sp
use statevarsmod,  only : invars,nd

implicit none

integer(i4),               intent(in)  :: year
integer(i4),               intent(in)  :: grid
real(sp),    dimension(:), intent(out) :: tmp        ! mean monthly temperature (degC)
real(sp),    dimension(:), intent(out) :: dtr        ! mean monthly diurnal temperature range (degC)
real(sp),    dimension(:), intent(out) :: pre        ! total monthly precipitation (mm)
real(sp),    dimension(:), intent(out) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
real(sp),    dimension(:), intent(out) :: cld        ! mean monthly cloud cover (fraction)
real(sp),    dimension(:), intent(out) :: wnd        ! mean monthly 10m windspeed (m s-1)
integer(i4), dimension(:), intent(out) :: ynd

integer(i4) :: start
integer(i4) :: end

!-------------------

start = (12 * year) + 1     ! Start month (i.e. Jan) begin after 12 months of the buffer year
end   = start + 11          ! End month (i.e. Dec)

start = start - 4           ! Include +/- 4 months buffer at beginning and end for gwgen calculations
end   = end + 4

!-------------------
! Copy monthly series from monvars into genvars for gwgen() input

tmp = invars(grid)%tmp(start:end)
dtr = invars(grid)%dtr(start:end)
pre = invars(grid)%pre(start:end)
wet = invars(grid)%wet(start:end)
cld = invars(grid)%cld(start:end)
wnd = invars(grid)%wnd(start:end)

ynd  = nd(start:end)

end subroutine copymonvars

!---------------------------------------------------------------------
!
! subroutine meanclim(rank,grid)
!
! ! Routine to calculate and store mean monthly climatology over the entire input time series
! ! Buffer months taken from start and last 4 months of mean climate
!
! use statevarsmod, only : nd,monvars,genvars_mean,calcyrs,cnt
! use weathergenmod, only : roundto
!
! implicit none
!
! integer(i4), intent(in) :: rank
! integer(i4), intent(in) :: grid
!
! real(sp), pointer, dimension(:) :: tmp        ! mean monthly temperature (degC)
! real(sp), pointer, dimension(:) :: dtr        ! mean monthly diurnal temperature range (degC)
! real(sp), pointer, dimension(:) :: pre        ! total monthly precipitation (mm)
! real(sp), pointer, dimension(:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
! real(sp), pointer, dimension(:) :: cld        ! mean monthly cloud cover (fraction)
! real(sp), pointer, dimension(:) :: wnd        ! mean monthly 10m windspeed (m s-1)
!
! integer :: gridcount
! integer :: aveyr
! integer :: yr
! integer :: i
! integer :: j
!
! !-------------------
!
! if (grid == 1) then
!   gridcount = cnt(1)
!   allocate(genvars_mean(gridcount))
! end if
!
! !-------------------
!
! tmp => genvars_mean(grid)%tmp
! dtr => genvars_mean(grid)%dtr
! pre => genvars_mean(grid)%pre
! wet => genvars_mean(grid)%wet
! cld => genvars_mean(grid)%cld
! wnd => genvars_mean(grid)%wnd
!
! !-------------------
!
! tmp = 0.
! dtr = 0.
! pre = 0.
! wet = 0.
! cld = 0.
! wnd = 0.
!
! aveyr = 20
!
! do yr = calcyrs-19, calcyrs
!   do i = 1, 20
!
!     j = i - 4
!
!     if (monvars(grid)%tmp(12*yr+j) < -273.15) monvars(grid)%tmp(12*yr+j) = 0.
!     if (monvars(grid)%dtr(12*yr+j) < 0.) monvars(grid)%dtr(12*yr+j) = 0.
!     if (monvars(grid)%pre(12*yr+j) < 0.) monvars(grid)%pre(12*yr+j) = 0.
!     if (monvars(grid)%wet(12*yr+j) < 0.) monvars(grid)%wet(12*yr+j) = 0.
!     if (monvars(grid)%cld(12*yr+j) < 0.) monvars(grid)%cld(12*yr+j) = 0.
!     if (monvars(grid)%wnd(12*yr+j) < 0.) monvars(grid)%wnd(12*yr+j) = 0.
!
!     tmp(i) = tmp(i) + monvars(grid)%tmp(12*yr+j) / real(aveyr)
!     dtr(i) = dtr(i) + monvars(grid)%dtr(12*yr+j) / real(aveyr)
!     pre(i) = pre(i) + monvars(grid)%pre(12*yr+j) / real(aveyr)
!     wet(i) = wet(i) + monvars(grid)%wet(12*yr+j) / real(aveyr)
!     cld(i) = cld(i) + monvars(grid)%cld(12*yr+j) / real(aveyr)
!     wnd(i) = wnd(i) + monvars(grid)%wnd(12*yr+j) / real(aveyr)
!
!   end do
! end do
!
! ! tmp = tmp / real(calcyrs)
! ! dtr = dtr / real(calcyrs)
! ! pre = pre / real(calcyrs)
! ! wet = wet / real(calcyrs)
! ! cld = cld / real(calcyrs)
! ! wnd = wnd / real(calcyrs)
! !
! ! do i = 5, 16
! !
! !   if (abs(tmp(i)) < 1.) tmp(i) = 0.
! !   if (abs(dtr(i)) < 1.) dtr(i) = 0.
! !   if (abs(pre(i)) < 1.) pre(i) = 0.
! !   if (abs(wet(i)) < 1.) wet(i) = 0.
! !   if (abs(cld(i)) < 1.) cld(i) = 0.
! !   if (abs(wnd(i)) < 1.) wnd(i) = 0.
! !
! ! end do
!
! tmp = roundto(tmp,1)
! dtr = roundto(dtr,1)
! pre = roundto(pre,1)
! wet = roundto(wet,1)
! cld = roundto(cld,1)
! wnd = roundto(wnd,2)
!
! genvars_mean(grid)%nd = [30,31,30,31, &
!                          31,28,31,30,31,30,31,31,30,31,30,31, &
!                          31,28,31,30]
!
! end subroutine meanclim
!
! !-----------------------------------------------------------------------

end module initmod
