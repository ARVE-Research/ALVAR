module drivermod

use parametersmod, only : i2,i4,sp,dp
use errormod,      only : ncstat,netcdf_err
use outputmod,     only : infompi
use netcdf
use mpi

implicit none

contains

!---------------------------------------------------------------------

subroutine initdate(info,job,rank)

use parametersmod, only : baseyr,ndaymonth
use metvarsmod   , only : startyr,endyr,calcyrs,nd,srt,cnt,cntt,p0,p1,tlen

implicit none

type(infompi), target    , intent(in) :: info
integer(i4), dimension(2), intent(in) :: job
integer                  , intent(in) :: rank

! Pointers for mpi info variables
character(100), pointer :: infile
character(100), pointer :: outfile
character(100), pointer :: timestring
integer(i4)   , pointer :: nproc
integer(i4)   , pointer :: t0
integer(i4)   , pointer :: nt

integer :: t1
integer :: yr
integer :: mon
integer :: i

!----------------------------------------------------
! Pointers to mpi derived type variables

infile     => info%infile
outfile    => info%outfile
timestring => info%timestring
nproc      => info%nproc
t0         => info%t0
nt         => info%nt

!----------------------------------------------------
! Get the buffered start and end of the time dimension of each validpixel

startyr = ((t0 - 1) / 12) + baseyr

calcyrs = nt / 12

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
  write(0,*)startyr,calcyrs
  write(0,*)cntt,nt
  write(0,*)t0,t1
  write(0,*)p0,p1
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

!-------------------
! Create start and count array from the mpi job info

srt = [job(1), t0]    ! Start of gridcell, start of time dimension
cnt = [job(2), nt]    ! Count of gridcell, count of time dimension


end subroutine initdate

!---------------------------------------------------------------------

subroutine initlonlat(info,job,rank)

! 1. Save current lon and lat of the gridcells into gridlon and gridlat
! 2. Decide if the current CPU process would require to print

use metvarsmod, only : lon,lat,indx,srt,cnt,gridlon,gridlat,lprint,gprint

implicit none

type(infompi), target    , intent(in) :: info
integer(i4), dimension(2), intent(in) :: job
integer                  , intent(in) :: rank

real(sp), pointer :: plon
real(sp), pointer :: plat

integer(i4) :: gridcount
integer(i4) :: ll
integer(i4) :: lon_loc
integer(i4) :: lat_loc
integer(i4), dimension(2) :: ll_loc

integer :: grid
integer :: i
integer :: j

!---------------------

plon => info%plon
plat => info%plat

gridcount = job(2)

allocate(gridlon(gridcount))
allocate(gridlat(gridcount))

!---------------------

do grid = 1, gridcount

  ll = srt(1) + (grid - 1)   ! Get the current index value

  ! Get the value of lon and lat from index dimension
  ll_loc = findloc(indx, ll)

  lon_loc = ll_loc(1)
  lat_loc = ll_loc(2)

  gridlon(grid) = lon(lon_loc)
  gridlat(grid) = lat(lat_loc)

end do

!---------------------

lprint = .FALSE.    ! logical variable to determine whether grid was sent to this core
gprint = 0          ! variable to store the grid index if lprint = .TRUE.

do grid = 1, gridcount

    if (gridlon(grid) == plon .AND. gridlat(grid) == plat) then

      lprint = .TRUE.
      gprint = grid

    end if

end do

! print *, lprint, plon, plat, gprint


end subroutine initlonlat

!---------------------------------------------------------------------

subroutine initmonvars()

use metvarsmod, only : cnt,cntt,monvars

implicit none

integer :: gridcount
integer :: moncount
integer :: i

!-------------------

gridcount = cnt(1)

!-------------------

allocate(monvars(gridcount))

do i = 1, gridcount

  allocate(monvars(i)%tmp(cntt))
  allocate(monvars(i)%dtr(cntt))
  allocate(monvars(i)%pre(cntt))
  allocate(monvars(i)%wet(cntt))
  allocate(monvars(i)%cld(cntt))
  allocate(monvars(i)%wnd(cntt))

end do


end subroutine initmonvars

!---------------------------------------------------------------------

subroutine copygenvars(year,grid)

use metvarsmod,    only : nd,monvars,genvars

implicit none

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid

real(sp), pointer, dimension(:) :: tmp        ! mean monthly temperature (degC)
real(sp), pointer, dimension(:) :: dtr        ! mean monthly diurnal temperature range (degC)
real(sp), pointer, dimension(:) :: pre        ! total monthly precipitation (mm)
real(sp), pointer, dimension(:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
real(sp), pointer, dimension(:) :: cld        ! mean monthly cloud cover (fraction)
real(sp), pointer, dimension(:) :: wnd        ! mean monthly 10m windspeed (m s-1)

integer(i4) :: start
integer(i4) :: end

!-------------------

tmp => genvars%tmp
dtr => genvars%dtr
pre => genvars%pre
wet => genvars%wet
cld => genvars%cld
wnd => genvars%wnd

!-------------------

start = (12 * year) + 1     ! Start month (i.e. Jan) begin after 12 months of the buffer year
end   = start + 11          ! End month (i.e. Dec)

start = start - 4           ! Include +/- 4 months buffer at beginning and end for gwgen calculations
end   = end + 4

!-------------------
! Copy monthly series from monvars into genvars for gwgen() input

tmp = monvars(grid)%tmp(start:end)
dtr = monvars(grid)%dtr(start:end)
pre = monvars(grid)%pre(start:end)
wet = monvars(grid)%wet(start:end)
cld = monvars(grid)%cld(start:end)
wnd = monvars(grid)%wnd(start:end)

genvars%nd  = nd(start:end)

end subroutine copygenvars

!---------------------------------------------------------------------
! Subroutine to allocate module variable 'dayvars' to either 365 or 366 days in a year

subroutine initdayvars(year,gridcount)

use metvarsmod,    only : nd,ndyear,dayvars

implicit none

integer(i4), intent(in) :: year
integer(i4), intent(in) :: gridcount

integer(i4) :: start
integer(i4) :: end

!-------------------

start = (12 * year) + 1     ! Start month (i.e. Jan) begin after 12 months of the buffer year
end   = start + 11          ! End month (i.e. Dec)

ndyear = sum(nd(start:end))

!-------------------
! Allocate number of days in year + 31 days (Jan of next year)
! Additional month needed for diurnal temp calculations

allocate(dayvars(gridcount,ndyear+31))


end subroutine initdayvars

!-----------------------------------------------------------------------

subroutine saveclonlat(grid)

! Save current lon and lat of the current gridcell into 'clon' and 'clat'

use metvarsmod,    only : lon,lat,indx,srt,cnt,clon,clat

implicit none

integer(i4), intent(in) :: grid

integer(i4) :: ll
integer(i4) :: lon_loc
integer(i4) :: lat_loc
integer(i4), dimension(2) :: ll_loc

!---------------------

ll = srt(1) + (grid - 1)   ! Get the current index value

!---------------------

! Get the value of lon and lat from index dimension
ll_loc = findloc(indx, ll)

lon_loc = ll_loc(1)
lat_loc = ll_loc(2)

clon = lon(lon_loc)
clat = lat(lat_loc)

end subroutine saveclonlat

!-----------------------------------------------------------------------

end module drivermod
