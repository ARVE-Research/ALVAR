module metvarsmod

use parametersmod, only : i2,i4,sp,dp
use orbitmod,      only : orbitpars

! This module contains the data structures that are read, calculated and stored by other subroutines
! Contains:
! --- full monthly array from a single netCDF file read (mon_metvars)
! --- full daily array of variables produced by gwgenmod (day_metvars)

!---------------------------------------------------------------------

type mon_metvars
  ! Derived datatype for monthly met variables input
  ! Dimension allocated from number of months

  real(sp), allocatable, dimension(:) :: tmp        ! mean monthly temperature (degC)
  real(sp), allocatable, dimension(:) :: dtr        ! mean monthly diurnal temperature range (degC)
  real(sp), allocatable, dimension(:) :: pre        ! total monthly precipitation (mm)
  real(sp), allocatable, dimension(:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
  real(sp), allocatable, dimension(:) :: cld        ! mean monthly cloud cover (fraction)
  real(sp), allocatable, dimension(:) :: wnd        ! mean monthly 10m windspeed (m s-1)

end type mon_metvars

type(mon_metvars), allocatable, dimension(:) :: monvars

!---------------------------------------------------------------------

type gwgen_vars
  ! Derived datatype storing 20 months of monthly varaibles (+/- 4 months buffer) for
  ! converting one year of data into daily values in gwgen()

  real(sp), dimension(20) :: tmp        ! mean monthly temperature (degC)
  real(sp), dimension(20) :: dtr        ! mean monthly diurnal temperature range (degC)
  real(sp), dimension(20) :: pre        ! total monthly precipitation (mm)
  real(sp), dimension(20) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
  real(sp), dimension(20) :: cld        ! mean monthly cloud cover (fraction)
  real(sp), dimension(20) :: wnd        ! mean monthly 10m windspeed (m s-1)

  integer , dimension(20) :: nd         ! number of days in month

end type gwgen_vars

type(gwgen_vars), target :: genvars

!---------------------------------------------------------------------

type day_metvars
  ! Derived datatype for the daily met variables output

  real(sp) :: prec    ! 24 hour total precipitation (mm)
  real(sp) :: tmin    ! 24 hour mean minimum temperature (degC)
  real(sp) :: tmax    ! 24 hour mean maximum temperature (degC)
  real(sp) :: cldf    ! 24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
  real(sp) :: wind    ! wind speed (m s-1)

  ! New added variables for PET calculations (Leo O Lai 16 Apr 2021)
  real(sp) :: tmean   ! 24 hour mean temperature (degC)
  real(sp) :: tday    ! mean daytime temperature (degC)
  real(sp) :: tnight  ! mean nighttime temperature (degC)
  real(sp) :: Pjj     ! precipitation equitability index for calculating PET
  real(sp) :: tdew    ! dew point temperature (degC)
  ! real(sp) :: tdew2
  real(sp) :: rhum    ! relative humidity (%)
  real(sp) :: dsol    ! solar declination angle (degree)
  real(sp) :: dayl    ! daylength (h)
  real(sp) :: srad    ! downwelling surface shortwave radiation (kJ m-2 d-1)
  real(sp) :: dpet    ! total potential evapotranspiration (mm)
  real(sp) :: daet    ! total actual evapotranspiration (mm)
  real(sp) :: alpha   ! ratio of daily AET/PET (fraction)
  real(sp) :: vpd     ! average daytime saturation vapor pressure deficit (Pa)

  real(dp), dimension(24) :: hprec    ! hourly precipitation (mm)

end type day_metvars

type(day_metvars), target, allocatable, dimension(:,:) :: dayvars ! Dimension allocate from number of days in year (365 or 366)

!---------------------------------------------------------------------

! Publicly shared variables that are constant across all gridcells
! Accessible to all modules and subroutines to avoid repetitive file read and calculations
! Used in initdate(), initmonvars() and initdayvars() subroutines
integer :: xlen                        ! length of dimension 'lat'
integer :: ylen                        ! length of dimension 'long'
integer :: ilen                        ! length of dimension 'index'
integer :: tlen                        ! length of dimension 'time'

! Allocatable arrays for longitude and latitude
real(dp),    allocatable, dimension(:)   :: lon
real(dp),    allocatable, dimension(:)   :: lat
real(dp),    allocatable, dimension(:)   :: time
integer(i4), allocatable, dimension(:,:) :: indx

! Module variable for storing lon and lat that user-specified for printing
real(dp) :: clon
real(dp) :: clat

logical :: lprint     ! TRUE if process recieved user-specified lon/lat and require printing
integer :: gprint     ! Save index value of the grid for printing to terminal

! Module variable for saving lon lat for each grid index
real(dp), allocatable, dimension(:) :: gridlon
real(dp), allocatable, dimension(:) :: gridlat

! Module variables for determining date and size of monthly data arrays
integer, allocatable, dimension(:) :: nd      ! number of days in months
integer :: startyr
integer :: endyr
integer :: calcyrs
integer :: ndyear          ! number of days in current year of calculation

! Module variables for reading in monthly series from netcdf file
integer, dimension(2) :: srt
integer, dimension(2) :: cnt
integer :: cntt
integer :: p0
integer :: p1


end module metvarsmod
