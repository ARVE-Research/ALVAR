module statevarsmod

use parametersmod, only : i2,i4,sp,dp
use orbitmod,      only : orbitpars

! This module contains the data structures that are read, calculated and stored by other subroutines

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
  real(sp) :: Ratm        ! Relative atmospheric pressure to sea-level (fraction)
  real(sp) :: Ratm30      ! Relative atmospheric pressure to 30m above sea-level (fraction)
  real(sp) :: Patm        ! Atmospheric pressure to sea-level (Pa)
  real(sp) :: Patm30      ! Atmospheric pressure to 30m above sea-level (Pa)

  real(sp) :: prec        ! 24 hour total precipitation (mm)
  real(sp) :: tmin        ! 24 hour mean minimum temperature (degC)
  real(sp) :: tmax        ! 24 hour mean maximum temperature (degC)
  real(sp) :: cldf        ! 24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
  real(sp) :: wind        ! wind speed (m s-1)

  ! New added variables (Leo O Lai 16 Apr 2021)
  real(sp) :: dayl           ! daylength (h)
  real(sp) :: dayl_n         ! daylength of the next day (h)
  integer(i4) :: sunrise     ! Current day sunrise hour (index of 24h hourly array)
  integer(i4) :: sunset      ! Current day sunset hour (index of 24h hourly array)
  integer(i4) :: sunrise_n   ! next day sunrise hour (index of 24h hourly array)
  integer(i4) :: dayhour     ! day time hours (h) --> from current day sunrise to sunset
  integer(i4) :: nighthour   ! night time hours (h) --> from current day sunset to next day sunrise

  real(sp) :: tmean       ! 24 hour mean temperature (degC)
  real(sp) :: tday        ! mean daytime temperature (degC)
  real(sp) :: tnight      ! mean nighttime temperature (degC)
  real(sp) :: Pjj         ! precipitation equitability index for calculating PET
  real(sp) :: tdew        ! dew point temperature (degC)
  real(sp) :: rhum        ! relative humidity (%)
  real(sp) :: dsol        ! solar declination angle (degree)
  real(sp) :: srad        ! downwelling surface shortwave radiation (kJ m-2 d-1)
  real(sp) :: srad_dir    ! direct beam downwelling shortwave raditaion (kJ m-2 d-1)
  real(sp) :: srad_dif    ! diffuse downwelling shortwave raditaion (kJ m-2 d-1)
  real(sp) :: lrad        ! upswelling surface longwave radiation (kJ m-2 d-1)
  real(sp) :: dpet        ! total potential evapotranspiration (mm)
  real(sp) :: daet        ! total actual evapotranspiration (mm)
  real(sp) :: alpha       ! ratio of daily AET/PET (fraction)
  real(sp) :: vpd         ! average daytime saturation vapor pressure deficit (Pa)

  ! Fire index variables
  real(sp) :: DF          ! Drought factor
  real(sp) :: KBDI        ! Keetch-Byram Drought Index (mm equivalent)
  real(sp) :: ForFireMk5  ! Forest fire danger index Mark 5 meter

  real(dp), dimension(24) :: hprec    ! hourly precipitation (mm)

end type day_metvars

type(day_metvars), target, allocatable, dimension(:,:) :: dayvars ! Dimension allocate from number of days in year (365 or 366)

!---------------------------------------------------------------------

! Number of soil layers
integer(i4), parameter :: nl = 6

type soildata

  logical :: validcell

  ! Derived datatype for soil layer variables
  real(sp), dimension(nl) :: sand       ! Sand content by mass (percent)
  real(sp), dimension(nl) :: clay       ! Clay content by mass (percent)
  real(sp), dimension(nl) :: cfvo       ! Course fragment content by volume (percent) --> Vcf variable in ARVE-DGVM
  real(sp), dimension(nl) :: OrgM       ! Organic matter content by mass (percent)

  real(sp), dimension(nl) :: Vsand      ! Sand content by volume (fraction)
  real(sp), dimension(nl) :: Vclay      ! Clay content by volume (fraction)
  real(sp), dimension(nl) :: Vsilt      ! Silt content by volume (fraction)
  real(sp), dimension(nl) :: VOrgM      ! Organic matter content by volume (fraction)
  real(sp), dimension(nl) :: rock       ! Course fragment content by mass (fraction)
  real(sp), dimension(nl) :: bulk       ! Soil bulk density (kg m-3)

  real(sp), dimension(nl) :: whc        ! Soil water holding capcity / available water content (mm)
  real(sp), dimension(nl) :: Ksat       ! Soil water saturated conductivity (mm s-1)
  real(sp), dimension(nl) :: Tsat       ! Soil water volumetric water content at saturation (fraction / m3 m-3)
  real(sp), dimension(nl) :: Tpor       ! Soil volumetric porosity (liquid minus ice content (or stones?)) (fraction)
  real(sp), dimension(nl) :: Ffrz       ! Fractional impermeable area as a function of soil ice content at a layer (fraction)
  real(sp), dimension(nl) :: Tfield     ! Soil water volumetric content at field capacity (Psi = -33 kPa)   (fraction / m3 m-3)
  real(sp), dimension(nl) :: Twilt      ! Soil water volumetric content at wilting point (Psi = -1500 kPa)   (fraction / m3 m-3)
  real(sp), dimension(nl) :: Psat       ! Soil water matric potential at saturation (mm)
  real(sp), dimension(nl) :: Bexp       ! Soil water B exponent used in the Brooks & Corey Pedotransfer Functions

  real(sp), dimension(nl) :: Csolid     ! Soil solids volumetric heat capcity at layer midpoint (J m-3 K-1)
  real(sp), dimension(nl) :: Ksolid     ! Soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
  ! real(sp), dimension(nl) :: Kdrysolid  !
  real(sp), dimension(nl) :: Kdry       ! Soil thermal conductivity of dry natural soil (W m-1 K-1)
  real(sp), dimension(nl) :: Kl         ! Soil thermal conductivity across layer boundary (W m-1 K-1)
  real(sp), dimension(nl) :: Cl         ! Instantaneous soil volumetric heat capacity at layer midpoint (J m-3 K-1)
  real(sp), dimension(nl) :: Kh         ! Snow thermal conductivity at layer midpoint (W m-1 K-1)

  real(sp), dimension(nl) :: Wliq       ! Soil liquid water content at layer midpoint (mm)
  real(sp), dimension(nl) :: Wice       ! Soil ice content at layer midpoint (mm)
  real(sp), dimension(nl) :: Tsoil      ! Soil temperature (K)
  real(sp), dimension(nl) :: Tsoiln     ! Soil temperature for precious timestep (K)
  real(sp), dimension(nl) :: Tliq       ! Soil volumetric water content (fraction / m3 m-3) / Fractional soil water content
  real(sp), dimension(nl) :: Tice       ! Soil volumetric ice content (fraction / m3 m-3) / Fraction soil ice content
  real(sp), dimension(nl) :: Psi        ! Soil water potential
  real(sp), dimension(nl) :: Psi_eq     ! Restriction for min of soil potential (mm)
  real(sp), dimension(nl) :: Ku         ! Soil water instantaneous (unsaturated) conductivity across layer boundary (mm s-1)

  real(sp) :: zw                        ! Mean water table depth (dimensionless)
  real(sp) :: fsat                      ! Saturated fraction / partial contributing area
  real(sp) :: Waquif_a                  ! Water in unconfined aquifer below soil (mm)
  real(sp) :: Waquif_t                  ! Water in aquifer within soil column (mm).
  real(sp) :: dTliq_b
  real(sp) :: raw_b                     ! Aerodynamic resistance for latent heat transfer btwn bare ground & atmosphere (s m-1)
  real(sp) :: rah_b                     ! Aerodynamic resistance for sensible heat transfer btwn bare ground & atmosphere (s m-1)

end type soildata

type(soildata), target, allocatable, dimension(:) :: soilvars   ! Dimension allocate to number of gridcells

!---------------------------------------------------------------------

type veg_vars

  ! Derived datatype for vegetation variables
  real(sp) :: gammastar         ! Photorespiratory compensation point (Pa)
  real(sp) :: gpp0              ! Gross primary productivity under non-water stressed condition (g C m-2 d-1) --> (Wang et al., 2017)
  real(sp) :: gpp               ! Gross primary productivity under actual condition (g C m-2 d-1) --> (Sitch et al., 2003; LPJ)
  real(sp) :: gpp_tot           ! Total daily gross primary productivity (g C d-1)
  real(sp) :: npp               ! Net primary productivity (g C m-2 d-1)
  real(sp) :: npp_tot           ! Total net primary productivity (g C d-1)
  real(sp) :: aresp             ! Autotrophic maintenence respiration (g C m-2 d-1)
  real(sp) :: rd                ! Daily leaf respiration (gC m-2 day-1) >> whole day include day + night
  real(sp) :: chi               ! Actual leaf internal / external CO2 partial pressure ratio (fraction)
  real(sp) :: chi0              ! Optimal leaf internal / external CO2 partial pressure ratio (fraction)
  real(sp) :: dgp               ! Optimal daily canopy conductance (mm m-2 s-1)
  real(sp) :: dgc               ! Actual daily canopy conductance (mm m-2 s-1
  real(sp) :: lai               ! LAI from data input (m2 m-2)
  real(sp) :: sla               ! Specific leaf area (m2 gC-1)

  integer(i4) :: biome          ! Biome classification from BIOME1 subroutine

  ! LPJ adapted vegetation variables
  logical  :: present           ! PFT present
  logical  :: estab             ! PFT establishment
  logical  :: survive           ! PFT survival
  real(sp) :: dwscal            ! Daily water stress factor (supply/demand ratio)
  real(sp) :: fpc_grid          ! Foilage projective cover over grid (fraction)
  real(sp) :: fpc_ind           ! Foliage projective cover of individual (fraction)
  real(sp) :: fpc_inc           ! Foliage projective cover increment (fraction)
  real(sp) :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
  real(sp) :: rm_ind            ! Root carbon mass of individual (gC m-2)
  real(sp) :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
  real(sp) :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
  real(sp) :: nind              ! PFT population
  real(sp) :: stemdiam          ! Tree stem diameter (m)
  real(sp) :: height            ! Tree height (m)
  real(sp) :: crownarea         ! Tree crownarea (m2)
  real(sp) :: lai_ind           ! Leaf area index of individual (m2 m-2)
  real(sp) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
  real(sp) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
  real(sp) :: litter_bg         ! Below ground litter pool (gC m-2)
  real(sp) :: turnover_ind      ! Total turnover of individual (gC m-2)

end type veg_vars

type(veg_vars), target, allocatable, dimension(:,:) :: vegvars   ! Dimension allocate to number of gridcells

!---------------------------------------------------------------------

type topo_vars
  ! Derived datatype for topographic variables

  real(sp) :: elev      ! Elevation (m)
  real(sp) :: slope     ! Slope (deg)
  real(sp) :: aspect    ! Aspect (deg)
  real(sp) :: cellarea  ! Area of gridcell (m2)
  real(sp) :: areafrac  ! Ground area fraction in gridcell (fraction)

  real(sp) :: sloperad      ! slope ratio

end type topo_vars

type(topo_vars), target, allocatable, dimension(:) :: topovars    ! Dimension allocate by number of gridcells

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

! Allocatable arrays for elevation (m)
real(sp), allocatable, dimension(:) :: elev

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


end module statevarsmod
