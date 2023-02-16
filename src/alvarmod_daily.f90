module alvarmod_daily

! Module containing the main ALVAR model

implicit none

real :: aaet

contains

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

subroutine alvar_daily(yr,grid,ndyear,calcyrs,day,i,spinup)

use parametersmod,    only : i4,sp,dp,Tfreeze
use statevarsmod,     only : sv,ns,nl
use pftparmod,        only : tree,evergreen,summergreen,raingreen,needle,boreal,c4
use randomdistmod,    only : randomstate
use radfluxmod,       only : radiativeflux
use latsensfluxmod,   only : latsensflux
use soilphysicsmod,   only : soilthermalprop,soilresistance,soiltemperature
use gppmod,           only : calcgpp_opt,calcgpp
use nppmod,           only : calcnpp
use waterbalancemod,  only : waterbalance
use soilhydrologymod, only : soilwaterflux
use snowhydrologymod, only : newsnowflux,snowdynamics
use lightmod,         only : light
use turnovermod,      only : turnover
use allocationmod,    only : allocation
use rootmod,          only : rootdist
use killplantmod,     only : killplant
use mortalitymod,     only : mortality
use averagemod,       only : longtermave

implicit none

integer(i4), intent(in) :: yr
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: ndyear
integer(i4), intent(in) :: calcyrs
integer(i4), intent(in) :: day
integer(i4), intent(in) :: i
logical,     intent(in) :: spinup

! Pointers for statevars
real(dp),          pointer :: lon
real(dp),          pointer :: lat
type(randomstate), pointer :: georndst

real(sp),    pointer, dimension(:) :: tmp
real(sp),    pointer, dimension(:) :: dtr
real(sp),    pointer, dimension(:) :: pre
real(sp),    pointer, dimension(:) :: wet
real(sp),    pointer, dimension(:) :: cld
real(sp),    pointer, dimension(:) :: wnd
integer(i4), pointer, dimension(:) :: nd

real(sp),    pointer               :: tmin         ! 24 hour mean minimum temperature (degC)
real(sp),    pointer               :: tmax         ! 24 hour mean maximum temperature (degC)
real(sp),    pointer               :: prec         ! 24 hour total precipitation (mm)
real(sp),    pointer               :: cldf         ! 24 hour mean cloud cover (fraction)
real(sp),    pointer               :: wind         ! 24 hour mean wind speed (m s-1)
real(sp),    pointer               :: dayl         ! Daylength (h)
integer(i4), pointer               :: sunrise      ! Current day sunrise hour (index of 24h hourly array)
integer(i4), pointer               :: sunset       ! Current day sunset hour (index of 24h hourly array)
integer(i4), pointer               :: sunrise_n    ! Next day sunrise hour (index of 24h hourly array)
integer(i4), pointer               :: dayhour      ! Day time hours (h) --> from current day sunrise to sunset
integer(i4), pointer               :: nighthour    ! Night time hours (h) --> from current day sunset to next day sunrise
real(sp),    pointer               :: tmean        ! 24 hour mean temperature (degC)
real(sp),    pointer               :: tday         ! Mean daytime temperature (degC)
real(sp),    pointer               :: tnight       ! Mean nighttime temperature (degC)
real(sp),    pointer               :: Ratm         ! Relative atmospheric pressure to sea-level (fraction)
real(sp),    pointer               :: Ratm30       ! Relative atmospheric pressure to 30m above sea-level (fraction)
real(sp),    pointer               :: Patm         ! Atmospheric pressure to sea-level (Pa)
real(sp),    pointer               :: Patm30       ! Atmospheric pressure to 30m above sea-level (Pa)
real(sp),    pointer               :: tdew         ! Dew point temperature (degC)
real(sp),    pointer               :: Pjj          ! Precipitation equitability index for calculating PET
real(sp),    pointer               :: rhum         ! Relative humidity (%)
real(sp),    pointer               :: dsol         ! Solar declination angle (degree)
real(sp),    pointer               :: srad         ! Downwelling surface shortwave radiation (kJ m-2 d-1)
real(sp),    pointer               :: srad_dir     ! Direct beam downwelling shortwave raditaion (kJ m-2 d-1)
real(sp),    pointer               :: srad_dif     ! Diffuse downwelling shortwave raditaion (kJ m-2 d-1)
real(sp),    pointer               :: lrad         ! Upswelling surface longwave radiation (kJ m-2 d-1)
real(sp),    pointer               :: dpet         ! Total potential evapotranspiration (mm)
real(sp),    pointer               :: daet         ! Total actual evapotranspiration (mm)
real(sp),    pointer               :: alpha        ! Ratio of daily AET/PET (fraction)
real(sp),    pointer               :: vpd          ! Average daytime saturation vapor pressure deficit (Pa)
real(sp),    pointer               :: dayprec
real(sp),    pointer               :: nightprec
real(dp),    pointer, dimension(:) :: hprec        ! Hourly precipitation (mm)
real(dp),    pointer, dimension(:) :: hprec_n      ! Next day hourly precipitation (mm)

logical,  pointer                 :: validcell
integer(i4), pointer              :: snl
real(sp), pointer, dimension(:)   :: dz              ! Snow/soil layer thickness (m)
real(sp), pointer, dimension(:)   :: dzmm            ! Snow/soil layer thickness (mm)
real(sp), pointer, dimension(:)   :: zpos            ! Midpoint z position (depth) of soil layer
real(sp), pointer, dimension(:)   :: zposmm          ! zpos in mm
real(sp), pointer, dimension(:)   :: zipos           ! Snow/soil layer interface z position (depth), positive downwards (dim = nl+1, including surface and column bottom interface)
real(sp), pointer, dimension(:)   :: ziposmm         ! zipos in mm
real(sp), pointer, dimension(:)   :: sand            ! Sand content by mass (percent)
real(sp), pointer, dimension(:)   :: clay            ! Clay content by mass (percent)
real(sp), pointer, dimension(:)   :: cfvo            ! Course fragment content by volume (percent) --> Vcf variable in ARVE-DGVM
real(sp), pointer, dimension(:)   :: OrgM            ! Organic matter content by mass (percent)
real(sp), pointer, dimension(:)   :: rock            ! Course fragment content by mass (fraction)
real(sp), pointer, dimension(:)   :: bulk            ! Soil bulk density (kg m-3)
real(sp), pointer, dimension(:)   :: Vsand           ! Sand content by volume (fraction)
real(sp), pointer, dimension(:)   :: Vclay           ! Clay content by volume (fraction)
real(sp), pointer, dimension(:)   :: Vsilt           ! Silt content by volume (fraction)
real(sp), pointer, dimension(:)   :: VOrgM           ! Organic matter content by volume (fraction)
real(sp), pointer, dimension(:)   :: awc             ! Soil water holding capcity / available water content (mm)
real(sp), pointer, dimension(:)   :: Ksat            ! Soil water saturated conductivity (mm s-1)
real(sp), pointer, dimension(:)   :: Tsat            ! Soil water volumetric water content at saturation (fraction / m3 m-3)
real(sp), pointer, dimension(:)   :: Tpor            ! Soil volumetric porosity (liquid minus ice content (or stones?)) (fraction)
real(sp), pointer, dimension(:)   :: Tfield          ! Soil water volumetric content at field capacity (Psi = -33 kPa)   (fraction / m3 m-3)
real(sp), pointer, dimension(:)   :: Twilt           ! Soil water volumetric content at wilting point (Psi = -1500 kPa)   (fraction / m3 m-3)
real(sp), pointer, dimension(:)   :: Psat            ! Soil water matric potential at saturation (mm)
real(sp), pointer, dimension(:)   :: Bexp            ! Soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(sp), pointer, dimension(:)   :: Csolid          ! Soil solids volumetric heat capcity at layer midpoint (J m-3 K-1)
real(sp), pointer, dimension(:)   :: Ksolid          ! Soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
real(sp), pointer, dimension(:)   :: Kdry            ! Soil thermal conductivity of dry natural soil (W m-1 K-1)
real(sp), pointer, dimension(:)   :: Kl              ! Soil thermal conductivity across layer boundary (W m-1 K-1)
real(sp), pointer, dimension(:)   :: Cl              ! Instantaneous soil volumetric heat capacity at layer midpoint (J m-3 K-1)
real(sp), pointer, dimension(:)   :: Kh              ! Snow thermal conductivity at layer midpoint (W m-1 K-1)
real(sp), pointer, dimension(:)   :: fact            ! Factor used in computing heat flux tridiagonal coefficients
real(sp), pointer, dimension(:)   :: Fhti            ! Heat flux across soil layer boundary (W m-2)
real(sp), pointer, dimension(:)   :: swf             ! Water status (fraction of whc/awc) --> for water routine adapted from LPJ-LMFire, variable "w" in waterbalanceLPJ.f90
real(sp), pointer, dimension(:)   :: aetsoil         ! Soil layer total actual evapotranspiration for all PFTs (mm)
real(sp), pointer, dimension(:)   :: Wliq            ! Soil liquid water content at layer midpoint (mm)
real(sp), pointer, dimension(:)   :: Wice            ! Soil ice content at layer midpoint (mm)
real(sp), pointer, dimension(:)   :: fice0           ! Layer ice fraction, previous timestep
real(sp), pointer, dimension(:)   :: Tsoil           ! Soil temperature (K)
real(sp), pointer, dimension(:)   :: Tsoiln          ! Soil temperature for precious timestep (K)
real(sp), pointer, dimension(:)   :: Tliq            ! Soil volumetric water content (fraction / m3 m-3) / Fractional soil water content
real(sp), pointer, dimension(:)   :: Tice            ! Soil volumetric ice content (fraction / m3 m-3) / Fraction soil ice content
real(sp), pointer, dimension(:)   :: Psi             ! Soil water potential
real(sp), pointer, dimension(:)   :: Psi_eq          ! Restriction for min of soil potential (mm)
real(sp), pointer, dimension(:)   :: Ku              ! Soil water instantaneous (unsaturated) conductivity across layer boundary (mm s-1)
real(sp), pointer                 :: zw              ! Mean water table depth (dimensionless)
real(sp), pointer                 :: fsat            ! Saturated fraction / partial contributing area
real(sp), pointer                 :: surf_runoff
real(sp), pointer                 :: surf_infl
real(sp), pointer                 :: Wliq_surf
real(sp), pointer                 :: Waquif_a        ! Water in unconfined aquifer below soil (mm)
real(sp), pointer                 :: Waquif_t        ! Water in aquifer within soil column (mm).
real(sp), pointer                 :: dTliq_b
real(sp), pointer                 :: raw_b           ! Aerodynamic resistance for latent heat transfer btwn bare ground & atmosphere (s m-1)
real(sp), pointer                 :: rah_b           ! Aerodynamic resistance for sensible heat transfer btwn bare ground & atmosphere (s m-1)
real(sp), pointer                 :: Wsno            ! Snow water equivalent of the snowpack (mm)
real(sp), pointer                 :: Wsno_old        ! Snow water equivalent of the snowpack (mm) from previous timestep
real(sp), pointer                 :: zsno            ! Snow total depth (m)
real(sp), pointer                 :: tausno          ! Snow age (non-dimensional)
real(sp), pointer                 :: tausno_old      ! Snow age of previous timestep
real(sp), pointer                 :: qsnomelt        ! Snow melt (kg m-2 s-1)
real(sp), pointer, dimension(:)   :: Ffrz            ! Fractional impermeable area as a function of soil ice content at a layer (fraction)
integer(i4), pointer, dimension(:)                 :: ithaw

real(sp), pointer                 :: toa_sw          ! Top of the atmosphere downwelling shortwave rad (kJ m-2 d-1)
real(sp), pointer                 :: Latm            ! Downwelling atmospheric longwave radiation (W m-2)
real(sp), pointer, dimension(:)   :: Lwg_v           ! Canopy net longwave rad flux for the ground (W m-2)
real(sp), pointer                 :: Lwg_b           ! Downwelling atmospheric longwave radiation into the surface (W m-2) bare surfaces
real(sp), pointer                 :: dLgdT_b         ! Derivative of the longwave radiation flux into the surface (W m-2) w.r.t. Temp
real(sp), pointer, dimension(:)   :: dLgdt_v         ! Derivative of the longwave radiation flux into the surface (W m-2) w.r.t. veg temp
real(sp), pointer, dimension(:)   :: Lwvd            ! Downward longwave below vegetation (W m-2)
real(sp), pointer, dimension(:)   :: Lwv             ! Net longwave rad flux for vegetation (W m-2)
real(sp), pointer, dimension(:)   :: dLvdT           ! Vegetation derivative of the longwave radiation flux into surface (W m-2) w,r,t Temp
real(sp), pointer                 :: Hg              ! Total sensible heat flux into the surface (W m-2)
real(sp), pointer                 :: dHgdT           ! Total derivative of the sensible heat flux into the surface (W m-2) w.r.t. temp
real(sp), pointer                 :: Hg_b            ! Sensible heat flux into the surface (W m-2) bare surfaces
real(sp), pointer                 :: dHgdt_b         ! Derivative of the sensible heat flux into the surface (W m-2) w.r.t temp (bare surf)
real(sp), pointer, dimension(:)   :: Hg_v            ! Sensible heat flux into the surface (W m-2) vegetated surfaces
real(sp), pointer, dimension(:)   :: dHgdt_v         ! Derivative of the sensible heat flux into surface (W m-2) w.r.t temp (vegetated surf)
real(sp), pointer                 :: Eg              ! Latent heat flux into the surface (W m-2)
real(sp), pointer                 :: dEgdT           ! Derivative of the latent heat flux into the surface (W m-2) with respect to Temperature
real(sp), pointer                 :: Eg_b            ! Latent heat flux into the surface (W m-2)(bare surfaces)
real(sp), pointer                 :: dEgdt_b         ! Derivative of the latent heat flux into surface (W m-2) w.r.t temp (bare surf)
real(sp), pointer, dimension(:)   :: Eg_v            ! Latent heat flux into the surface (W m-2)(vegetated surfaces)
real(sp), pointer, dimension(:)   :: dEgdt_v         ! Derivative of the latent heat flux into the surface (W m-2) w.r.t temp (veg surf)
real(sp), pointer                 :: Beta_soil       ! CLM 4.0 Eqn 5.68 represent molecular diffusion processes in evaporation from soil
real(sp), pointer                 :: dqsath_dT       ! Saturated water vapour specific humidity (CLM 4 Eqn 5.143) derivative w.r.t grnd temp.
real(sp), pointer                 :: dqgdTg          ! Derivative of the specific humidity of the soil surface (CLM 4.0 Eqn 5.74)
real(sp), pointer                 :: qg              ! Ground specific humidity
real(sp), pointer                 :: qs              ! Canopy specific humidity
real(sp), pointer                 :: lam             ! Latent heat of phase change at current time step (J kg-1)
real(sp), pointer                 :: hs              ! Net energy flux into the surface (W m-2)
real(sp), pointer                 :: dhsdT           ! Derivative of hs with respect to temperature
real(sp), pointer                 :: fsnow           ! Fraction of the gridcell covered by snow (fraction)
real(sp), pointer                 :: mo_runoff_tot   ! Monthly total overland runoff (mm)
real(sp), pointer                 :: qseva_tot       ! Timestep soil evaporation (mm)
real(sp), pointer                 :: qsubl_tot       ! Timestep soil sublimation (mm)
real(sp), pointer                 :: qsdew_tot       ! Timestep soil dew  (mm)
real(sp), pointer                 :: qfrost_tot      ! Timestep soil frost (mm)
real(sp), pointer                 :: qseva           ! Soil evaporation flux (mm s-1)
real(sp), pointer                 :: qsubl           ! Soil sublimation flux (mm s-1)
real(sp), pointer                 :: qsdew           ! Soil surface dew flux (mm s-1)
real(sp), pointer                 :: qfrost          ! Soil surface frost flux (mm s-1)
real(sp), pointer                 :: qliq            ! Liquid precipitation (kg m-2 sec-1)
real(sp), pointer                 :: qover           ! Liquid water surface runoff (kg m-2 sec-1)
real(sp), pointer                 :: qsno            ! Snowfall (kg m-2 sec -1)
real(sp), pointer, dimension(:)   :: qgrnd_l         ! Total rate of liquid precip reaching the ground under canopy(mm s-1)
real(sp), pointer, dimension(:)   :: qgrnd_s         ! Total rate of solid precip reaching the ground under canopy(mm s-1)

real(sp), pointer, dimension(:)   :: gammastar       ! Photorespiratory compensation point (Pa)
real(sp), pointer, dimension(:)   :: lambda          ! Actual leaf internal / external CO2 partial pressure ratio (fraction)
real(sp), pointer, dimension(:)   :: chi0            ! Optimal leaf internal / external CO2 partial pressure ratio (fraction)
real(sp), pointer, dimension(:)   :: gpp0            ! Gross primary productivity under non-water stressed condition (g C m-2 d-1)
real(sp), pointer, dimension(:)   :: gpp0_tot        ! Total optimal daily gross primary productivity (g C d-1)
real(sp), pointer, dimension(:)   :: gpp             ! Gross primary productivity under actual condition (g C m-2 d-1) --> (Sitch et al., 2003; LPJ)
real(sp), pointer, dimension(:)   :: gpp_tot         ! Total daily gross primary productivity (g C d-1)
real(sp), pointer, dimension(:)   :: npp             ! Net primary productivity (g C m-2 d-1)
real(sp), pointer, dimension(:)   :: npp_tot         ! Total net primary productivity (g C d-1)
real(sp), pointer, dimension(:)   :: aresp           ! Autotrophic maintenence respiration (g C m-2 d-1)
real(sp), pointer, dimension(:)   :: dgp             ! Optimal daily canopy conductance (mm m-2 s-1)
real(sp), pointer, dimension(:)   :: rd              ! Daily leaf respiration (gC m-2 day-1) >> whole day include day + night
real(sp), pointer, dimension(:)   :: dgc             ! Actual daily canopy conductance (mm m-2 s-1)
real(sp), pointer, dimension(:,:) :: dwscal          ! Daily water stress factor (supply/demand ratio)
real(sp), pointer, dimension(:,:) :: aet             ! Actual evapotranspiration per PFT (mm d-1)
real(sp), pointer, dimension(:)   :: dphen           ! Phenology status of summergreen (proportion of leaf-on) (fraction)
real(sp), pointer, dimension(:)   :: dphen_t         ! Temperature based phenology of summergreen (proportion of leaf-on) (fraction)
real(sp), pointer, dimension(:)   :: dphen_w         ! Water based phenology of summergreen (proportion of leaf-on) (fraction)

logical,  pointer, dimension(:)   :: present         ! PFT present
logical,  pointer, dimension(:)   :: estab           ! PFT establishment
logical,  pointer, dimension(:)   :: survive         ! PFT survival
logical,  pointer, dimension(:)   :: leafon          ! Leaf phenology on/off for photosynthesis
real(sp), pointer, dimension(:)   :: lai             ! LAI from data input (m2 m-2)
real(sp), pointer, dimension(:)   :: sla             ! Specific leaf area (m2 gC-1)
real(sp), pointer, dimension(:)   :: fpc_grid        ! Foilage projective cover over grid (fraction)
real(sp), pointer, dimension(:)   :: fpc_ind         ! Foliage projective cover of individual (fraction)
real(sp), pointer, dimension(:)   :: fpc_inc         ! Foliage projective cover increment (fraction)
real(sp), pointer, dimension(:)   :: lm_ind          ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:)   :: rm_ind          ! Root carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:)   :: sm_ind          ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:)   :: hm_ind          ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:)   :: clabile         ! Short-term carbon labile storage for growth respiration (gC m-2)
real(sp), pointer, dimension(:)   :: creserve        ! Long-term carbon reserve (gC m-2)
real(sp), pointer, dimension(:)   :: nind            ! PFT population
real(sp), pointer, dimension(:)   :: stemdiam        ! Tree stem diameter (m)
real(sp), pointer, dimension(:)   :: height          ! Tree height (m)
real(sp), pointer, dimension(:)   :: crownarea       ! Tree crownarea (m2)
real(sp), pointer, dimension(:)   :: lai_ind         ! Leaf area index of individual (m2 m-2)
real(sp), pointer, dimension(:)   :: litter_ag_fast  ! Fast above ground litter pool (gC m-2)
real(sp), pointer, dimension(:)   :: litter_ag_slow  ! Slow above ground litter pool (gC m-2)
real(sp), pointer, dimension(:)   :: litter_bg       ! Below ground litter pool (gC m-2)
real(sp), pointer, dimension(:)   :: turnover_ind    ! Total turnover of individual (gC m-2)
real(sp), pointer, dimension(:)   :: bm_inc          ! Total biomass increment (yearly at the moment) (gC m-2 yr-1)
real(sp), pointer, dimension(:)   :: leafondays      ! Number of days since leaf phenology is on
real(sp), pointer, dimension(:)   :: leafoffdays     ! Number of days since leaf phenology is off
real(sp), pointer, dimension(:)   :: meanfpc
real(sp), pointer, dimension(:)   :: abm_inc
real(sp), pointer, dimension(:)   :: Wcan            ! Amount of water in the canopy (mm)
real(sp), pointer, dimension(:)   :: Wcanmax         ! Maximum quantity of water the canopy can hold (mm)
real(sp), pointer, dimension(:)   :: can_evap        ! Canopy evaporation (mm s-1)
real(sp), pointer, dimension(:)   :: xs_can_drip     ! Excess canopy water (mm)
real(sp), pointer, dimension(:)   :: frac            ! Fraction of sapwood and heartwood that is below ground (in coarse roots)(gC)
real(sp), pointer, dimension(:)   :: rootdepth       ! Root depth (m)
real(sp), pointer, dimension(:,:) :: rootfracl       ! Root fraction in that soil layer

real(sp), pointer                 :: elev            ! Elevation (m)
real(sp), pointer                 :: slope           ! Slope (deg)
real(sp), pointer                 :: aspect          ! Aspect (deg)
real(sp), pointer                 :: cellarea        ! Area of gridcell (m2)
real(sp), pointer                 :: areafrac        ! Ground area fraction in gridcell (fraction)

!-------------------------

real(sp), dimension(12) :: mtmp
real(sp), dimension(12) :: mpre
real(sp) :: srad_W
real(sp) :: srad_dir_W
real(sp) :: srad_dif_W

integer(i4) :: d
integer(i4) :: pft

!-----------------------------------------------------------

validcell   => sv(grid)%validcell
lon         => sv(grid)%lon
lat         => sv(grid)%lat

tmp => sv(grid)%monvars%tmp
dtr => sv(grid)%monvars%dtr
pre => sv(grid)%monvars%pre
wet => sv(grid)%monvars%wet
cld => sv(grid)%monvars%cld
wnd => sv(grid)%monvars%wnd
nd  => sv(grid)%monvars%nd

snl         => sv(grid)%soilvars%snl
dz          => sv(grid)%soilvars%dz
dzmm        => sv(grid)%soilvars%dzmm
zpos        => sv(grid)%soilvars%zpos
zposmm      => sv(grid)%soilvars%zposmm
zipos       => sv(grid)%soilvars%zipos
ziposmm     => sv(grid)%soilvars%ziposmm
cfvo        => sv(grid)%soilvars%cfvo
Vsand       => sv(grid)%soilvars%Vsand
VOrgM       => sv(grid)%soilvars%VOrgM
awc         => sv(grid)%soilvars%awc
Ksat        => sv(grid)%soilvars%Ksat
Tsat        => sv(grid)%soilvars%Tsat
Tpor        => sv(grid)%soilvars%Tpor
Ffrz        => sv(grid)%soilvars%Ffrz
Tfield      => sv(grid)%soilvars%Tfield
Twilt       => sv(grid)%soilvars%Twilt
Psat        => sv(grid)%soilvars%Psat
Bexp        => sv(grid)%soilvars%Bexp
Csolid      => sv(grid)%soilvars%Csolid
Ksolid      => sv(grid)%soilvars%Ksolid
Kdry        => sv(grid)%soilvars%Kdry
Kl          => sv(grid)%soilvars%Kl
Cl          => sv(grid)%soilvars%Cl
Kh          => sv(grid)%soilvars%Kh
fact        => sv(grid)%soilvars%fact
Fhti        => sv(grid)%soilvars%Fhti
swf         => sv(grid)%soilvars%swf
aetsoil     => sv(grid)%soilvars%aetsoil
Wliq        => sv(grid)%soilvars%Wliq
Wice        => sv(grid)%soilvars%Wice
fice0       => sv(grid)%soilvars%fice0
Tsoil       => sv(grid)%soilvars%Tsoil
Tsoiln      => sv(grid)%soilvars%Tsoiln
Tliq        => sv(grid)%soilvars%Tliq
Tice        => sv(grid)%soilvars%Tice
Psi         => sv(grid)%soilvars%Psi
Psi_eq      => sv(grid)%soilvars%Psi_eq
Ku          => sv(grid)%soilvars%Ku
zw          => sv(grid)%soilvars%zw
fsat        => sv(grid)%soilvars%fsat
surf_runoff => sv(grid)%soilvars%surf_runoff
surf_infl   => sv(grid)%soilvars%surf_infl
Wliq_surf   => sv(grid)%soilvars%Wliq_surf
dTliq_b     => sv(grid)%soilvars%dTliq_b
raw_b       => sv(grid)%soilvars%raw_b
rah_b       => sv(grid)%soilvars%rah_b

Wsno          => sv(grid)%soilvars%Wsno
Wsno_old      => sv(grid)%soilvars%Wsno_old
zsno          => sv(grid)%soilvars%zsno
tausno        => sv(grid)%soilvars%tausno
tausno_old    => sv(grid)%soilvars%tausno_old
qsnomelt      => sv(grid)%soilvars%qsnomelt
ithaw         => sv(grid)%soilvars%ithaw

toa_sw        => sv(grid)%surfvars%toa_sw
Latm          => sv(grid)%surfvars%Latm
Lwg_v         => sv(grid)%surfvars%Lwg_v
Lwg_b         => sv(grid)%surfvars%Lwg_b
dLgdT_b       => sv(grid)%surfvars%dLgdT_b
dLgdt_v       => sv(grid)%surfvars%dLgdt_v
Lwvd          => sv(grid)%surfvars%Lwvd
Lwv           => sv(grid)%surfvars%Lwv
dLvdT         => sv(grid)%surfvars%dLvdT
Hg            => sv(grid)%surfvars%Hg
dHgdT         => sv(grid)%surfvars%dHgdT
Hg_b          => sv(grid)%surfvars%Hg_b
dHgdt_b       => sv(grid)%surfvars%dHgdt_b
Hg_v          => sv(grid)%surfvars%Hg_v
dHgdt_v       => sv(grid)%surfvars%dHgdt_v
Eg            => sv(grid)%surfvars%Eg
dEgdT         => sv(grid)%surfvars%dEgdT
Eg_b          => sv(grid)%surfvars%Eg
dEgdt_b       => sv(grid)%surfvars%dEgdt_b
Eg_v          => sv(grid)%surfvars%Eg_v
dEgdt_v       => sv(grid)%surfvars%dEgdt_v
Beta_soil     => sv(grid)%surfvars%Beta_soil
dqsath_dT     => sv(grid)%surfvars%dqsath_dT
dqgdTg        => sv(grid)%surfvars%dqgdTg
qg            => sv(grid)%surfvars%qg
qs            => sv(grid)%surfvars%qs
lam           => sv(grid)%surfvars%lam
hs            => sv(grid)%surfvars%hs
dhsdT         => sv(grid)%surfvars%dhsdT

fsnow         => sv(grid)%surfvars%fsnow
mo_runoff_tot => sv(grid)%surfvars%mo_runoff_tot
qseva_tot     => sv(grid)%surfvars%qseva_tot
qsubl_tot     => sv(grid)%surfvars%qsubl_tot
qsdew_tot     => sv(grid)%surfvars%qsdew_tot
qfrost_tot    => sv(grid)%surfvars%qfrost_tot
qseva         => sv(grid)%surfvars%qseva
qsubl         => sv(grid)%surfvars%qsubl
qsdew         => sv(grid)%surfvars%qsdew
qfrost        => sv(grid)%surfvars%qfrost
qliq          => sv(grid)%surfvars%qliq
qover         => sv(grid)%surfvars%qover
qsno          => sv(grid)%surfvars%qsno
qgrnd_l       => sv(grid)%surfvars%qgrnd_l
qgrnd_s       => sv(grid)%surfvars%qgrnd_s

tmin      => sv(grid)%dayvars%tmin(day)
tmax      => sv(grid)%dayvars%tmax(day)
prec      => sv(grid)%dayvars%prec(day)
cldf      => sv(grid)%dayvars%cldf(day)
wind      => sv(grid)%dayvars%wind(day)
dayl      => sv(grid)%dayvars%dayl(day)
dpet      => sv(grid)%dayvars%dpet(day)
daet      => sv(grid)%dayvars%daet(day)
sunrise   => sv(grid)%dayvars%sunrise(day)
sunset    => sv(grid)%dayvars%sunset(day)
sunrise_n => sv(grid)%dayvars%sunrise_n(day)
dayhour   => sv(grid)%dayvars%dayhour(day)
nighthour => sv(grid)%dayvars%nighthour(day)
tmean     => sv(grid)%dayvars%tmean(day)
tday      => sv(grid)%dayvars%tday(day)
tnight    => sv(grid)%dayvars%tnight(day)
tdew      => sv(grid)%dayvars%tdew(day)
Ratm      => sv(grid)%dayvars%Ratm(day)
Patm      => sv(grid)%dayvars%Patm(day)
Patm30    => sv(grid)%dayvars%Patm30(day)
Pjj       => sv(grid)%dayvars%Pjj(day)
rhum      => sv(grid)%dayvars%rhum(day)
srad      => sv(grid)%dayvars%srad(day)
srad_dir  => sv(grid)%dayvars%srad_dir(day)
srad_dif  => sv(grid)%dayvars%srad_dif(day)
dsol      => sv(grid)%dayvars%dsol(day)
lrad      => sv(grid)%dayvars%lrad(day)
hprec     => sv(grid)%dayvars%hprec(day,:)
hprec_n   => sv(grid)%dayvars%hprec(day+1,:)
dayprec   => sv(grid)%dayvars%dayprec(day)
nightprec => sv(grid)%dayvars%nightprec(day)

gammastar => sv(grid)%gppvars%gammastar(day,:)
lambda    => sv(grid)%gppvars%lambda(day,:)
chi0      => sv(grid)%gppvars%chi0(day,:)
gpp0      => sv(grid)%gppvars%gpp0(day,:)
gpp0_tot  => sv(grid)%gppvars%gpp_tot(day,:)
gpp       => sv(grid)%gppvars%gpp(day,:)
gpp_tot   => sv(grid)%gppvars%gpp_tot(day,:)
npp       => sv(grid)%gppvars%npp(day,:)
npp_tot   => sv(grid)%gppvars%npp_tot(day,:)
aresp     => sv(grid)%gppvars%aresp(day,:)
dgp       => sv(grid)%gppvars%dgp(day,:)
rd        => sv(grid)%gppvars%rd(day,:)
dgc       => sv(grid)%gppvars%dgc(day,:)
dwscal    => sv(grid)%gppvars%dwscal(:,:)
aet       => sv(grid)%gppvars%aet(:,:)
dphen     => sv(grid)%gppvars%dphen(day,:)
dphen_t   => sv(grid)%gppvars%dphen_t(day,:)
dphen_w   => sv(grid)%gppvars%dphen_w(day,:)

present        => sv(grid)%vegvars%present
estab          => sv(grid)%vegvars%estab
survive        => sv(grid)%vegvars%survive
leafon         => sv(grid)%vegvars%leafon
lai            => sv(grid)%vegvars%lai
sla            => sv(grid)%vegvars%sla
fpc_grid       => sv(grid)%vegvars%fpc_grid
fpc_ind        => sv(grid)%vegvars%fpc_ind
fpc_inc        => sv(grid)%vegvars%fpc_inc
lm_ind         => sv(grid)%vegvars%lm_ind
rm_ind         => sv(grid)%vegvars%rm_ind
sm_ind         => sv(grid)%vegvars%sm_ind
hm_ind         => sv(grid)%vegvars%hm_ind
clabile        => sv(grid)%vegvars%clabile
creserve       => sv(grid)%vegvars%creserve
nind           => sv(grid)%vegvars%nind
stemdiam       => sv(grid)%vegvars%stemdiam
height         => sv(grid)%vegvars%height
crownarea      => sv(grid)%vegvars%crownarea
lai_ind        => sv(grid)%vegvars%lai_ind
litter_ag_fast => sv(grid)%vegvars%litter_ag_fast
litter_ag_slow => sv(grid)%vegvars%litter_ag_slow
litter_bg      => sv(grid)%vegvars%litter_bg
turnover_ind   => sv(grid)%vegvars%turnover_ind
bm_inc         => sv(grid)%vegvars%bm_inc
leafondays     => sv(grid)%vegvars%leafondays
leafoffdays    => sv(grid)%vegvars%leafoffdays
meanfpc        => sv(grid)%vegvars%meanfpc
abm_inc        => sv(grid)%vegvars%abm_inc
Wcan           => sv(grid)%vegvars%Wcan
Wcanmax        => sv(grid)%vegvars%Wcanmax
can_evap       => sv(grid)%vegvars%can_evap
xs_can_drip    => sv(grid)%vegvars%xs_can_drip
frac           => sv(grid)%vegvars%frac
rootdepth      => sv(grid)%vegvars%rootdepth
rootfracl      => sv(grid)%vegvars%rootfracl

elev     => sv(grid)%topovars%elev
slope    => sv(grid)%topovars%slope
aspect   => sv(grid)%topovars%aspect
cellarea => sv(grid)%topovars%cellarea
areafrac => sv(grid)%topovars%areafrac

!-----------------------------------------------------------

mtmp = tmp(5:16)
mpre = pre(5:16)

srad_W      = (srad*1000.)     / (dayl*3600.)
srad_dir_W  = (srad_dir*1000.) / (dayl*3600.)
srad_dif_W  = (srad_dif*1000.) / (dayl*3600.)

swf = max(min(1.,Wliq(1:6)/awc),0.)

!-------------------------

call newsnowflux(i,dayl,present,snl,Wsno,zsno,tausno,tmean,prec,dz,zpos,zipos,Wliq,Wice,Tsoil,Tsoiln, &
                 fpc_grid,xs_can_drip,fsnow,qliq,qsno,qgrnd_l,qgrnd_s)

if (i == 1) call light(tree,present,fpc_grid,fpc_ind,fpc_inc,lm_ind,rm_ind,sm_ind,hm_ind,nind,sla,crownarea,lai_ind,&
                       litter_ag_fast,litter_ag_slow,litter_bg,turnover_ind,meanfpc)

if (i == 1) call calcgpp_opt(day,dayl,srad,srad_dir,tmean,zsno,tree,raingreen,summergreen,c4,present,dwscal,fpc_grid,&
                             cellarea,areafrac,leafon,leafondays,leafoffdays,dphen,dphen_t,dphen_w,dgp,rd,gpp0,gpp0_tot)

if (i == 1) call waterbalance(dpet,awc,swf,dgp,dphen,present,rm_ind,fpc_grid,rootfracl,daet,aet(day,:),aetsoil,dwscal(day,:),dgc)

call soilwaterflux(i,snl,validcell,dayl,daet,sunrise,sunset,sunrise_n,dayhour,nighthour,hprec,hprec_n,dz,dzmm,zpos,zposmm,&
                   zipos,ziposmm,dayprec,nightprec,awc,Ksat,Tsat,Tpor,Ffrz,Psat,Bexp,aetsoil,aet(day,:),fpc_grid,rootfracl,&
                   Wliq,Wice,Tsoil,Tsoiln,Tliq,Tice,Psi,Psi_eq,Ku,zw,fsat,surf_runoff,surf_infl,Wliq_surf,dTliq_b)

if (-snl > 0) call snowdynamics(i,dayl,snl,Wsno,Wsno_old,zsno,tausno,tausno_old,dz,zpos,zipos, &
                                Wliq,Wice,Tsoil,Tsoiln,fice0,ithaw,fsnow)

if (i == 1) call calcgpp(day,dayl,srad,srad_dir,tmean,c4,present,fpc_grid,cellarea,areafrac,dgc,lambda,rd,gpp,gpp_tot)

call calcnpp(day,i,dayl,tday,tnight,Tsoil(1),tree,present,gpp,dphen,leafon,leafondays,leafoffdays,lm_ind,rm_ind,sm_ind,hm_ind, &
             nind,cellarea,areafrac,clabile,creserve,aresp,npp,npp_tot,bm_inc,abm_inc)

call soilthermalprop(validcell,snl,dz,zpos,zipos,cfvo,Vsand,VOrgM,Tsat,Psat,Bexp,Tsoil,Tliq,Tice,Wliq,Wice,Wsno,&
                     Csolid,Ksolid,Kdry,Kl,Cl,Kh)

! call radiativeflux(day,i,snl,dayhour,nighthour,dayl,lat,cldf,tdew,dsol,srad,srad_dir,srad_dif)

call soilresistance(wind,raw_b,rah_b)

call latsensflux(i,snl,dpet,tday,tnight,tdew,cldf,wind,rhum,Patm,Patm30,srad_W,fsnow,zsno,raw_b,rah_b,Tfield,&
                 Tsoil,Tsoiln,Psat,Tsat,dz,Bexp,Wliq,Wice,&
                 Lwg_b,dLgdt_b,Hg,dHgdT,Hg_b,dHgdt_b,Eg,dEgdT,Eg_b,dEgdt_b,Beta_soil,dqsath_dT,dqgdTg,qg,qs,lam)

                 if (i == 1) then          ! Day time

                   if (dayl > 0.) then

                     hs = (srad)*1000./(dayl*3600.) - Lwg_b - Hg_b - lam * Eg_b  !5.129
                     dhsdT =  - dLgdT_b - dHgdT_b - lam * dEgdT_b

                   else

                     hs = 0. - Lwg_b - Hg_b - lam * Eg_b  !5.129
                     dhsdT =  - dLgdT_b - dHgdT_b - lam * dEgdT_b

                   end if

                 else if (i == 2) then     ! Night time

                   hs = 0. - Lwg_b - Hg_b - lam * Eg_b  !5.129
                   dhsdT =  - dLgdT_b - dHgdT_b - lam * dEgdT_b

                   ! if (dayl == 23.) then
                   !
                   !   hs = -(srad)*1000./(dayl*3600.)/24 - Lwg_b - Hg_b - lam * Eg_b  !5.129
                   !   dhsdT =  - dLgdT_b - dHgdT_b - lam * dEgdT_b
                   !
                   ! end if

                 end if

call soiltemperature(i,snl,Patm,dayl,dayhour,nighthour,tday,tnight,srad_dir,lrad,dz,dzmm,zpos,zipos,&
                     Wliq,Wice,Tsat,Bexp,Psat,Tpor,Tliq,Tice,fice0,Wsno,zsno,qsnomelt,fsnow,Kl,Cl,Kh, &
                     fact,Fhti,ithaw,raw_b,rah_b,hs,dhsdT,Tsoil,Tsoiln)

if (i == 2) call turnover(present,nind,lm_ind,rm_ind,sm_ind,hm_ind,litter_ag_fast,litter_ag_slow,litter_bg,turnover_ind)

if (i == 2) call allocation(day,tree,evergreen,present,dwscal,dphen,npp,bm_inc,fpc_grid,fpc_ind,fpc_inc,nind,sla,stemdiam,&
                            height,crownarea,lai_ind,lm_ind,rm_ind,sm_ind,hm_ind,clabile,litter_ag_fast,litter_ag_slow,litter_bg)

if (i == 2) call rootdist(tree,present,zipos,lm_ind,rm_ind,sm_ind,hm_ind,nind,crownarea,frac,rootdepth,rootfracl)

if (i == 2) call light(tree,present,fpc_grid,fpc_ind,fpc_inc,lm_ind,rm_ind,sm_ind,hm_ind,nind,sla,crownarea,lai_ind,&
                       litter_ag_fast,litter_ag_slow,litter_bg,turnover_ind,meanfpc)

if (yr > calcyrs - 30 .and. i == 2) call longtermave(30.,grid,day,ndyear)

if (day == ndyear .and. i == 2) then

  call killplant(abm_inc,present,lm_ind,rm_ind,sm_ind,hm_ind,nind,litter_ag_fast,litter_ag_slow,litter_bg)

  call mortality(mtmp,sv(grid)%dayvars%tmean,abm_inc,sla,present,lm_ind,rm_ind,sm_ind,hm_ind,nind,litter_ag_fast,litter_ag_slow,&
                 litter_bg,turnover_ind)

end if



!-----------------------------------------------------------


! print *,yr,day,i,dpet,daet,sum(aetsoil),sum(gpp0),sum(gpp),sum(npp),bm_inc,sum(height)
! print *, yr,day,i,prec,sum(aet(:,4)),sum(aet(:,5)),sum(aet(:,6)),sum(aet(:,7)),rootdepth
! print *, yr,day,i,prec,aet(day,4),aet(day,5),aet(day,6),aet(day,7),dpet,daet,tmean
! print *, yr,day,i,dayl,hs,(srad)*1000./(dayl*3600.),-Lwg_b,-Hg_b,-lam*Eg_b,rah_b
! print *, yr,day,i,dayl,hs,sum(meanfpc),(Tliq+Tice), sum(sv(grid)%gppvars%gpp), snl,zsno
! print *, yr,day,present,fpc_grid,(Tliq(1:6)+Tice(1:6))/Tsat(1:6), sum(sv(grid)%gppvars%gpp), prec, dpet,daet
! print *, yr,day,present,fpc_grid,(Tliq(1:6))/Tsat(1:6), sum(sv(grid)%gppvars%gpp),prec,dpet,daet,zsno,lon
! print *, yr,day,present,fpc_grid, sum(meanfpc),sum(sv(grid)%gppvars%gpp),sum(sv(grid)%gppvars%npp),prec,dpet,daet,zsno,lon
! print *, yr,day,elev,validcell,lon
! print *, yr,day,i,dayl, tday,tnight, hs,sum(meanfpc),Tsoil-Tfreeze, snl,zsno,sum(sv(grid)%gppvars%gpp)
! print *, yr,day,i,dayl, tday,tnight, Wice, snl,zsno,sum(sv(grid)%gppvars%gpp)
! print *, yr,day,i,hs,dayl,Tsoil(ns:nl)-Tfreeze,snl,zsno, sum(Wice), sum(meanfpc),lon
! pft = 4
! print *, yr,day,i,meanfpc,dphen_w(pft),dphen_t(pft),npp(pft),leafon(pft),leafondays(pft),leafoffdays(pft),tmean,zsno,&
!         creserve(pft),clabile(pft)
! print *, yr,day,i,rootfracl(3,:)
! print *, yr,day,i,rootfracl(pft,:)

! print *, '   '
! print *, dz, Tliq

! if (yr > 1020 .and. i == 2) print *, yr-1020,day,gpp,npp,sum(fpc_grid(1:7)),sum(fpc_grid(8:9)),elev


! if (day == ndyear .and. i == 2) then
!
!   print *, yr,sum(sv(grid)%gppvars%gpp(:,5)), &
!               sum(sv(grid)%gppvars%gpp(:,6)), &
!               sum(sv(grid)%gppvars%gpp(:,7)), &
!               sum(sv(grid)%gppvars%gpp(:,8)), &
!               sum(sv(grid)%gppvars%gpp(:,:)), &
!               sum(sv(grid)%gppvars%npp(:,5)), &
!               sum(sv(grid)%gppvars%npp(:,6)), &
!               sum(sv(grid)%gppvars%npp(:,7)), &
!               sum(sv(grid)%gppvars%npp(:,8)), &
!               sum(sv(grid)%gppvars%npp(:,:)), &
!               sum(fpc_grid(1:7)),sum(fpc_grid(8:9))
!
! end if
!
! print *, yr, day, sv(grid)%longave%gpp,sv(grid)%longave%npp















!-------------------------

! mtmp = tmp(5:16)
! mpre = pre(5:16)

!-------------------------




end subroutine alvar_daily

!-------------------------------------------------------

end module alvarmod_daily
