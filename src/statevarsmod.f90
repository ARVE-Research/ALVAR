module statevarsmod

use parametersmod, only : i2,i4,sp,dp
use pftparmod,     only : npft
use orbitmod,      only : orbitpars
use randomdistmod, only : randomstate

! This module contains the data structures that are read, calculated and stored by other subroutines

!---------------------------------------------------------------------

type in_metvars
  ! Derived datatype for monthly met variables input
  ! Dimension allocated from number of months

  real(sp), allocatable, dimension(:) :: tmp        ! Mean monthly temperature (degC)
  real(sp), allocatable, dimension(:) :: dtr        ! Mean monthly diurnal temperature range (degC)
  real(sp), allocatable, dimension(:) :: pre        ! Total monthly precipitation (mm)
  real(sp), allocatable, dimension(:) :: wet        ! Number of days in the month with precipitation > 0.1 mm (days)
  real(sp), allocatable, dimension(:) :: cld        ! Mean monthly cloud cover (fraction)
  real(sp), allocatable, dimension(:) :: wnd        ! Mean monthly 10m windspeed (m s-1)

end type in_metvars

type(in_metvars), allocatable, dimension(:) :: invars

!---------------------------------------------------------------------

type mon_metvars
  ! Derived datatype storing 20 months of monthly varaibles (+/- 4 months buffer) for
  ! converting one year of data into daily values in gwgen()

  real(sp),    dimension(20) :: tmp        ! Mean monthly temperature (degC)
  real(sp),    dimension(20) :: dtr        ! Mean monthly diurnal temperature range (degC)
  real(sp),    dimension(20) :: pre        ! Total monthly precipitation (mm)
  real(sp),    dimension(20) :: wet        ! Number of days in the month with precipitation > 0.1 mm (days)
  real(sp),    dimension(20) :: cld        ! Mean monthly cloud cover (fraction)
  real(sp),    dimension(20) :: wnd        ! Mean monthly 10m windspeed (m s-1)
  integer(i4), dimension(20) :: nd         ! Number of days in month

end type mon_metvars

! type(gwgen_vars), target :: genvars
! type(gwgen_vars), target, allocatable, dimension(:) :: genvars_mean

!---------------------------------------------------------------------

type day_metvars

  ! Derived datatype for the daily met variables output
  ! Dimension allocated as 366 + 31 = 397 days (number of days per year + days in Jan)
  real(sp), dimension(397) :: Ratm            ! Relative atmospheric pressure to sea-level (fraction)
  real(sp), dimension(397) :: Ratm30          ! Relative atmospheric pressure to 30m above sea-level (fraction)
  real(sp), dimension(397) :: Patm            ! Atmospheric pressure to sea-level (Pa)
  real(sp), dimension(397) :: Patm30          ! Atmospheric pressure to 30m above sea-level (Pa)97
  real(sp), dimension(397) :: tmin            ! 24 hour mean minimum temperature (degC)
  real(sp), dimension(397) :: tmax            ! 24 hour mean maximum temperature (degC)
  real(sp), dimension(397) :: prec            ! 24 hour total precipitation (mm)
  real(sp), dimension(397) :: cldf            ! 24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
  real(sp), dimension(397) :: wind            ! Wind speed (m s-1)

  ! New added variables (Leo O Lai 16 Apr 2021)
  real(sp),    dimension(397) :: dayl         ! Daylength (h)
  real(sp),    dimension(397) :: dayl_n       ! Daylength of the next day (h)
  integer(i4), dimension(397) :: sunrise      ! Current day sunrise hour (index of 24h hourly array)
  integer(i4), dimension(397) :: sunset       ! Current day sunset hour (index of 24h hourly array)
  integer(i4), dimension(397) :: sunrise_n    ! Next day sunrise hour (index of 24h hourly array)
  integer(i4), dimension(397) :: dayhour      ! Day time hours (h) --> from current day sunrise to sunset
  integer(i4), dimension(397) :: nighthour    ! Night time hours (h) --> from current day sunset to next day sunrise

  real(sp),    dimension(397) :: tmean        ! 24 hour mean temperature (degC)
  real(sp),    dimension(397) :: tday         ! Mean daytime temperature (degC)
  real(sp),    dimension(397) :: tnight       ! Mean nighttime temperature (degC)
  real(sp),    dimension(397) :: Pjj          ! Precipitation equitability index for calculating PET
  real(sp),    dimension(397) :: tdew         ! Dew point temperature (degC)
  real(sp),    dimension(397) :: rhum         ! Relative humidity (%)
  real(sp),    dimension(397) :: dsol         ! Solar declination angle (degree)
  real(sp),    dimension(397) :: srad         ! Downwelling surface shortwave radiation (kJ m-2 d-1)
  real(sp),    dimension(397) :: srad_dir     ! Direct beam downwelling shortwave raditaion (kJ m-2 d-1)
  real(sp),    dimension(397) :: srad_dif     ! Diffuse downwelling shortwave raditaion (kJ m-2 d-1)
  real(sp),    dimension(397) :: lrad         ! Upswelling surface longwave radiation (kJ m-2 d-1)
  real(sp),    dimension(397) :: dpet         ! Total potential evapotranspiration (mm)
  real(sp),    dimension(397) :: daet         ! Total actual evapotranspiration (mm)
  real(sp),    dimension(397) :: alpha        ! Ratio of daily AET/PET (fraction)
  real(sp),    dimension(397) :: vpd          ! Average daytime saturation vapor pressure deficit (Pa)

  real(sp),    dimension(397)    :: dayprec
  real(sp),    dimension(397)    :: nightprec
  real(dp),    dimension(397,24) :: hprec     ! Hourly precipitation (mm)

  ! Fire index variables
  real(sp),    dimension(397) :: DF           ! Drought factor
  real(sp),    dimension(397) :: KBDI         ! Keetch-Byram Drought Index (mm equivalent)
  real(sp),    dimension(397) :: ForFireMk5   ! Forest fire danger index Mark 5 meter

end type day_metvars

!---------------------------------------------------------------------

! Number of soil layers
integer(i4), parameter :: nl        = 6             ! Number of soil layers
integer(i4), parameter :: snomax    = 5             ! Maximum number of snow layers
integer(i4), parameter :: ns        = 1 - snomax    ! Index of top snow layer
real(sp),    parameter :: soildepth = 2.
real(sp),    parameter :: zsno_max  = 4.            ! maximum snow depth (m)

type soil_vars

  ! Derived datatype for soil and snow layer variables
  integer(i4)                  :: snl        ! Index value of top snow layer (negative is more layers)
  integer(i4)                  :: gnl        ! Index value of lowest soil layer
  real(sp), dimension(ns:nl)   :: dz         ! Snow/soil layer thickness (m)
  real(sp), dimension(ns:nl)   :: dzmm       ! Snow/soil layer thickness (mm)
  real(sp), dimension(ns:nl)   :: zpos       ! Midpoint z position (depth) of soil layer
  real(sp), dimension(ns:nl)   :: zposmm     ! zpos in mm
  real(sp), dimension(ns:nl+1) :: zipos      ! Snow/soil layer interface z position (depth), positive downwards (dim = nl+1, including surface and column bottom interface)
  real(sp), dimension(ns:nl+1) :: ziposmm    ! zipos in mm

  real(sp), dimension(nl)      :: sand       ! Sand content by mass (percent)
  real(sp), dimension(nl)      :: clay       ! Clay content by mass (percent)
  real(sp), dimension(nl)      :: cfvo       ! Course fragment content by volume (percent) --> Vcf variable in ARVE-DGVM
  real(sp), dimension(nl)      :: OrgM       ! Organic matter content by mass (percent)

  real(sp), dimension(nl)      :: Vsand      ! Sand content by volume (fraction)
  real(sp), dimension(nl)      :: Vclay      ! Clay content by volume (fraction)
  real(sp), dimension(nl)      :: Vsilt      ! Silt content by volume (fraction)
  real(sp), dimension(nl)      :: VOrgM      ! Organic matter content by volume (fraction)
  real(sp), dimension(nl)      :: rock       ! Course fragment content by mass (fraction)
  real(sp), dimension(nl)      :: bulk       ! Soil bulk density (kg m-3)

  real(sp), dimension(nl)      :: awc        ! Soil water holding capcity / available water content (mm)
  real(sp), dimension(ns:nl)   :: Ksat       ! Soil water saturated conductivity (mm s-1)
  real(sp), dimension(ns:nl)   :: Tsat       ! Soil water volumetric water content at saturation (fraction / m3 m-3)
  real(sp), dimension(ns:nl)   :: Tpor       ! Soil volumetric porosity (liquid minus ice content (or stones?)) (fraction)
  real(sp), dimension(nl)      :: Tfield     ! Soil water volumetric content at field capacity (Psi = -33 kPa)   (fraction / m3 m-3)
  real(sp), dimension(nl)      :: Twilt      ! Soil water volumetric content at wilting point (Psi = -1500 kPa)   (fraction / m3 m-3)
  real(sp), dimension(ns:nl)   :: Psat       ! Soil water matric potential at saturation (mm)
  real(sp), dimension(ns:nl)   :: Bexp       ! Soil water B exponent used in the Brooks & Corey Pedotransfer Functions

  real(sp), dimension(nl)      :: Csolid     ! Soil solids volumetric heat capcity at layer midpoint (J m-3 K-1)
  real(sp), dimension(nl)      :: Ksolid     ! Soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
  real(sp), dimension(nl)      :: Kdrysolid  !
  real(sp), dimension(nl)      :: Kdry       ! Soil thermal conductivity of dry natural soil (W m-1 K-1)
  real(sp), dimension(ns:nl)   :: Kl         ! Soil thermal conductivity across layer boundary (W m-1 K-1)
  real(sp), dimension(ns:nl)   :: Cl         ! Instantaneous soil volumetric heat capacity at layer midpoint (J m-3 K-1)
  real(sp), dimension(ns:nl)   :: Kh         ! Snow thermal conductivity at layer midpoint (W m-1 K-1)
  real(sp), dimension(ns:nl)   :: fact       ! Factor used in computing heat flux tridiagonal coefficients
  real(sp), dimension(ns:nl)   :: Fhti       ! Heat flux across soil layer boundary (W m-2)
  integer(i4), dimension(ns:nl) :: ithaw

  real(sp), dimension(nl)      :: swf        ! Water status (fraction of whc/awc) --> for water routine adapted from LPJ-LMFire, variable "w" in waterbalanceLPJ.f90
  real(sp), dimension(nl)      :: aetsoil    ! Soil layer total actual evapotranspiration for all PFTs (mm)
  real(sp), dimension(ns:nl)   :: fice0      ! Layer ice fraction, previous timestep
  real(sp), dimension(ns:nl)   :: Wliq       ! Soil liquid water content at layer midpoint (mm)
  real(sp), dimension(ns:nl)   :: Wice       ! Soil ice content at layer midpoint (mm)
  real(sp), dimension(ns:nl)   :: Tsoil      ! Soil temperature (K)
  real(sp), dimension(ns:nl)   :: Tsoiln     ! Soil temperature for precious timestep (K)
  real(sp), dimension(ns:nl)   :: Tliq       ! Soil volumetric water content (fraction / m3 m-3) / Fractional soil water content
  real(sp), dimension(ns:nl)   :: Tice       ! Soil volumetric ice content (fraction / m3 m-3) / Fraction soil ice content
  real(sp), dimension(ns:nl)   :: Psi        ! Soil water potential
  real(sp), dimension(ns:nl)   :: Psi_eq     ! Restriction for min of soil potential (mm)
  real(sp), dimension(ns:nl)   :: Ku         ! Soil water instantaneous (unsaturated) conductivity across layer boundary (mm s-1)

  ! Snow variables
  real(sp)                     :: Wsno       ! Snow water equivalent of the snowpack (mm)
  real(sp)                     :: Wsno_old   ! Snow water equivalent of the snowpack (mm) from previous timestep
  real(sp)                     :: zsno       ! Snow total depth (m)
  real(sp)                     :: tausno     ! Snow age (non-dimensional)
  real(sp)                     :: tausno_old ! Snow age of previous timestep
  real(sp), dimension(nl)      :: Ffrz       ! Fractional impermeable area as a function of soil ice content at a layer (fraction)
  real(sp)                     :: qsnomelt   ! Snow melt (kg m-2 s-1)

  real(sp) :: zw                        ! Mean water table depth (dimensionless)
  real(sp) :: fsat                      ! Saturated fraction / partial contributing area
  real(sp) :: Wliq_surf
  real(sp) :: surf_runoff               ! Surface water runoff (mm)
  real(sp) :: surf_infl                 ! Surface water infiltration (mm)
  real(sp) :: Waquif_a                  ! Water in unconfined aquifer below soil (mm)
  real(sp) :: Waquif_t                  ! Water in aquifer within soil column (mm).
  real(sp) :: dTliq_b
  real(sp) :: raw_b                     ! Aerodynamic resistance for latent heat transfer btwn bare ground & atmosphere (s m-1)
  real(sp) :: rah_b                     ! Aerodynamic resistance for sensible heat transfer btwn bare ground & atmosphere (s m-1)

end type soil_vars

!---------------------------------------------------------------------

integer(i4), parameter :: band = 2

type surf_vars

  real(sp) :: toa_sw                            ! Top of the atmosphere downwelling shortwave rad (kJ m-2 d-1)
  real(sp) :: sw
  real(sp) :: sw_dir
  real(sp) :: sw_dif

  !Longwave
  real(sp) :: Latm                              ! Downwelling atmospheric longwave radiation (W m-2)
  real(sp), dimension(npft) :: Lwg_v            ! Canopy net longwave rad flux for the ground (W m-2)
  real(sp) :: Lwg_b                             ! Downwelling atmospheric longwave radiation into the surface (W m-2) bare surfaces
  real(sp) :: dLgdT_b                           ! Derivative of the longwave radiation flux into the surface (W m-2) w.r.t. Temp
  real(sp), dimension(npft) :: dLgdt_v          ! Derivative of the longwave radiation flux into the surface (W m-2) w.r.t. veg temp
  real(sp), dimension(npft) :: Lwvd             ! Downward longwave below vegetation (W m-2)
  real(sp), dimension(npft) :: Lwv              ! Net longwave rad flux for vegetation (W m-2)
  real(sp), dimension(npft) :: dLvdT            ! Vegetation derivative of the longwave radiation flux into surface (W m-2) w,r,t Temp

  !Sensible heat fluxes
  real(sp) :: Hg                                ! Total sensible heat flux into the surface (W m-2)
  real(sp) :: dHgdT                             ! Total derivative of the sensible heat flux into the surface (W m-2) w.r.t. temp
  real(sp) :: Hg_b                              ! Sensible heat flux into the surface (W m-2) bare surfaces
  real(sp) :: dHgdt_b                           ! Derivative of the sensible heat flux into the surface (W m-2) w.r.t temp (bare surf)
  real(sp), dimension(npft) :: Hg_v             ! Sensible heat flux into the surface (W m-2) vegetated surfaces
  real(sp), dimension(npft) :: dHgdt_v          ! Derivative of the sensible heat flux into surface (W m-2) w.r.t temp (vegetated surf)

  !Latent heat fluxese
  real(sp) :: Eg                                ! Latent heat flux into the surface (W m-2)
  real(sp) :: dEgdT                             ! Derivative of the latent heat flux into the surface (W m-2) with respect to Temperature
  real(sp) :: Eg_b                              ! Latent heat flux into the surface (W m-2)(bare surfaces)
  real(sp) :: dEgdt_b                           ! Derivative of the latent heat flux into surface (W m-2) w.r.t temp (bare surf)
  real(sp), dimension(npft) :: Eg_v             ! Latent heat flux into the surface (W m-2)(vegetated surfaces)
  real(sp), dimension(npft) :: dEgdt_v          ! Derivative of the latent heat flux into the surface (W m-2) w.r.t temp (veg surf)
  real(sp) :: Beta_soil                         ! CLM 4.0 Eqn 5.68 represent molecular diffusion processes in evaporation from soil
  real(sp) :: dqsath_dT                         ! Saturated water vapour specific humidity (CLM 4 Eqn 5.143) derivative w.r.t grnd temp.
  real(sp) :: dqgdTg                            ! Derivative of the specific humidity of the soil surface (CLM 4.0 Eqn 5.74)
  real(sp) :: qg                                ! Ground specific humidity
  real(sp) :: qs                                ! Canopy specific humidity

  !Total flux
  real(sp) :: lam                               ! Latent heat of phase change at current time step (J kg-1)
  real(sp) :: hs                                ! Net energy flux into the surface (W m-2)
  real(sp) :: dhsdT                             ! Derivative of hs with respect to temperature

  ! Soil fluxes
  real(sp) :: fsnow                             ! Fraction of the gridcell covered by snow (fraction)
  real(sp) :: mo_runoff_tot                     ! Monthly total overland runoff (mm)
  real(sp) :: qseva_tot                         ! Timestep soil evaporation (mm)
  real(sp) :: qsubl_tot                         ! Timestep soil sublimation (mm)
  real(sp) :: qsdew_tot                         ! Timestep soil dew  (mm)
  real(sp) :: qfrost_tot                        ! Timestep soil frost (mm)

  real(sp) :: qseva                             ! Soil evaporation flux (mm s-1)
  real(sp) :: qsubl                             ! Soil sublimation flux (mm s-1)
  real(sp) :: qsdew                             ! Soil surface dew flux (mm s-1)
  real(sp) :: qfrost                            ! Soil surface frost flux (mm s-1)
  real(sp) :: qliq                              ! Liquid precipitation (kg m-2 sec-1)
  real(sp) :: qover                             ! Liquid water surface runoff (kg m-2 sec-1)
  real(sp) :: qsno                              ! Snowfall (kg m-2 sec -1)
  real(sp), dimension(npft) :: qgrnd_l          ! Total rate of liquid precip reaching the ground under canopy(mm s-1)
  real(sp), dimension(npft) :: qgrnd_s          ! Total rate of solid precip reaching the ground under canopy(mm s-1)

end type surf_vars

!---------------------------------------------------------------------

type gpp_vars

  ! Derived datatype for daily productivity variables (GPP and NPP)
  real(sp), dimension(366,npft) :: gammastar      ! Photorespiratory compensation point (Pa)
  real(sp), dimension(366,npft) :: lambda         ! Actual leaf internal / external CO2 partial pressure ratio (fraction)
  real(sp), dimension(366,npft) :: chi0           ! Optimal leaf internal / external CO2 partial pressure ratio (fraction)
  real(sp), dimension(366,npft) :: gpp0           ! Gross primary productivity under non-water stressed condition (g C m-2 d-1)
  real(sp), dimension(366,npft) :: gpp0_tot       ! Total optimal daily gross primary productivity (g C d-1)
  real(sp), dimension(366,npft) :: gpp            ! Gross primary productivity under actual condition (g C m-2 d-1) --> (Sitch et al., 2003; LPJ)
  real(sp), dimension(366,npft) :: gpp_tot        ! Total daily gross primary productivity (g C d-1)
  real(sp), dimension(366,npft) :: npp            ! Net primary productivity (g C m-2 d-1)
  real(sp), dimension(366,npft) :: npp_tot        ! Total net primary productivity (g C d-1)
  real(sp), dimension(366,npft) :: aresp          ! Autotrophic maintenence respiration (g C m-2 d-1)
  real(sp), dimension(366,npft) :: dgp            ! Optimal canopy conductance (mm m-2 s-1)
  real(sp), dimension(366,npft) :: rd             ! Daily leaf respiration (gC m-2 day-1) >> whole day include day + night
  real(sp), dimension(366,npft) :: dgc            ! Actual canopy conductance (mm m-2 s-1)
  real(sp), dimension(366,npft) :: dwscal         ! Daily water stress factor (supply/demand ratio)
  real(sp), dimension(366,npft) :: aet            ! Actual evapotranspiration per PFT (mm d-1)

  real(sp), dimension(366,npft) :: dphen          ! Phenology status of summergreen (proportion of leaf-on) (fraction)
  real(sp), dimension(366,npft) :: dphen_t        ! Temperature based phenology of summergreen (proportion of leaf-on) (fraction)
  real(sp), dimension(366,npft) :: dphen_w        ! Water based phenology of summergreen (proportion of leaf-on) (fraction)

end type gpp_vars

!---------------------------------------------------------------------

type veg_vars

  ! Derived datatype for vegetation variables adapted from LPJ-LMFire
  logical,  dimension(npft)    :: present          ! PFT present
  logical,  dimension(npft)    :: estab            ! PFT establishment
  logical,  dimension(npft)    :: survive          ! PFT survival
  logical,  dimension(npft)    :: leafon           ! Leaf phenology on/off for photosynthesis
  real(sp), dimension(npft)    :: lai              ! LAI from data input (m2 m-2)
  real(sp), dimension(npft)    :: sla              ! Specific leaf area (m2 gC-1)
  real(sp), dimension(npft)    :: fpc_grid         ! Foilage projective cover over grid (fraction)
  real(sp), dimension(npft)    :: fpc_ind          ! Foliage projective cover of individual (fraction)
  real(sp), dimension(npft)    :: fpc_inc          ! Foliage projective cover increment (fraction)
  real(sp), dimension(npft)    :: lm_ind           ! Leaf carbon mass of individual (gC m-2)
  real(sp), dimension(npft)    :: rm_ind           ! Root carbon mass of individual (gC m-2)
  real(sp), dimension(npft)    :: sm_ind           ! Sapwood carbon mass of individual (gC m-2)
  real(sp), dimension(npft)    :: hm_ind           ! Heartwood carbon mass of individual (gC m-2)
  real(sp), dimension(npft)    :: clabile          ! Short-term carbon labile storage for growth respiration (gC m-2)
  real(sp), dimension(npft)    :: creserve         ! Long-term carbon reserve (gC m-2)
  real(sp), dimension(npft)    :: nind             ! PFT population
  real(sp), dimension(npft)    :: stemdiam         ! Tree stem diameter (m)
  real(sp), dimension(npft)    :: height           ! Tree height (m)
  real(sp), dimension(npft)    :: crownarea        ! Tree crownarea (m2)
  real(sp), dimension(npft)    :: lai_ind          ! Leaf area index of individual (m2 m-2)
  real(sp), dimension(npft)    :: litter_ag_fast   ! Fast above ground litter pool (gC m-2)
  real(sp), dimension(npft)    :: litter_ag_slow   ! Slow above ground litter pool (gC m-2)
  real(sp), dimension(npft)    :: litter_bg        ! Below ground litter pool (gC m-2)
  real(sp), dimension(npft)    :: turnover_ind     ! Total turnover of individual (gC m-2)
  real(sp), dimension(npft)    :: bm_inc           ! Total biomass increment (yearly at the moment) (gC m-2 yr-1)
  real(sp), dimension(npft)    :: leafondays       ! Number of days since leaf phenology is on
  real(sp), dimension(npft)    :: leafoffdays      ! Number of days since leaf phenology is off
  real(sp), dimension(npft)    :: meanfpc
  real(sp), dimension(npft)    :: abm_inc

  real(sp), dimension(npft)    :: Wcan             ! Amount of water in the canopy (mm)
  real(sp), dimension(npft)    :: Wcanmax          ! Maximum quantity of water the canopy can hold (mm)
  real(sp), dimension(npft)    :: can_evap         ! Canopy evaporation (mm s-1)
  real(sp), dimension(npft)    :: xs_can_drip      ! Excess canopy water (mm)

  real(sp), dimension(npft)    :: frac             ! Fraction of sapwood and heartwood that is below ground (in coarse roots)(gC)
  real(sp), dimension(npft)    :: rootdepth        ! Root depth (m)
  real(sp), dimension(npft,nl) :: rootfracl        ! Root fraction in that soil layer

  integer(i4) :: biome                             ! Biome classification from BIOME1 subroutine

end type veg_vars

!---------------------------------------------------------------------

type topo_vars

  ! Derived datatype for topographic variables
  real(sp) :: elev          ! Elevation (m)
  real(sp) :: slope         ! Slope (deg)
  real(sp) :: aspect        ! Aspect (deg)
  real(sp) :: cellarea      ! Area of gridcell (m2)
  real(sp) :: areafrac      ! Ground area fraction in gridcell (fraction)
  real(sp) :: sloperad      ! Ratio of shortwave radiation on slope to flat surface (sloperad/flatrad)

end type topo_vars

!---------------------------------------------------------------------

type long_ave

  ! Derived type to save long term average of statevars for output
  real(sp) :: daet
  real(sp) :: gpp
  real(sp) :: npp
  real(sp) :: lm_ind
  real(sp) :: rm_ind
  real(sp) :: sm_ind
  real(sp) :: hm_ind
  real(sp) :: height
  real(sp) :: soilw
  real(sp) :: Tsoil
  real(sp) :: zsno

end type long_ave

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

type statevars

  type(mon_metvars) :: monvars
  type(day_metvars) :: dayvars
  type(soil_vars)   :: soilvars
  type(surf_vars)   :: surfvars
  type(gpp_vars)    :: gppvars
  type(veg_vars)    :: vegvars
  type(topo_vars)   :: topovars
  type(long_ave)    :: longave

  real(dp)          :: lon
  real(dp)          :: lat
  type(randomstate) :: georndst

  logical           :: validcell

end type statevars

type(statevars), target, allocatable, dimension(:) :: sv

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

! Publicly shared variables that are constant across all gridcells
! Accessible to all modules and subroutines to avoid repetitive file read and calculations
integer :: xlen                        ! length of dimension 'lat'
integer :: ylen                        ! length of dimension 'long'
integer :: ilen                        ! length of dimension 'index'
integer :: tlen                        ! length of dimension 'time'

! Allocatable arrays for longitude and latitude
real(dp),    allocatable, dimension(:)   :: inlon
real(dp),    allocatable, dimension(:)   :: inlat
real(dp),    allocatable, dimension(:)   :: time
integer(i4), allocatable, dimension(:,:) :: indx

! Module variables for determining date and size of monthly data arrays
integer(i4), allocatable, dimension(:) :: nd      ! number of days in each months
integer(i4) :: startyr
integer(i4) :: endyr
integer(i4) :: calcyrs

end module statevarsmod
