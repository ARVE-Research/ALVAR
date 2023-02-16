module alvarmod

! Module containing the main ALVAR model

implicit none

contains

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

subroutine alvar_annual(yr,grid,ndyear,calcyrs,spinup)

use parametersmod,    only : i4,sp,dp
use statevarsmod,     only : sv,invars
use randomdistmod,    only : randomstate,initrndstate
use orbitmod,         only : orbitpars,orbit
use soilstatemod,     only : initsoil
use gwgenmodnew,      only : gwgen_new
use diurnaltempmod,   only : diurnaltemp,humidity
use radiationmod,     only : elev_Ratm,calcPjj,radpet,calctdew,calcVPD
use hourlyprecmod,    only : hourlyprec
use establishmentmod, only : sapling,bioclim,establishment
use phenologymod,     only : summerphenology
use initmod,          only : copymonvars

implicit none

integer(i4), intent(in) :: yr
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: ndyear
integer(i4), intent(in) :: calcyrs
logical,     intent(in) :: spinup

! Pointers for statevars
real(dp),          pointer :: lon
real(dp),          pointer :: lat
real(sp),          pointer :: elev
type(randomstate), pointer :: georndst

real(sp),    pointer, dimension(:) :: tmp
real(sp),    pointer, dimension(:) :: dtr
real(sp),    pointer, dimension(:) :: pre
real(sp),    pointer, dimension(:) :: wet
real(sp),    pointer, dimension(:) :: cld
real(sp),    pointer, dimension(:) :: wnd
integer(i4), pointer, dimension(:) :: nd

real(sp),    pointer, dimension(:)   :: tmin         ! 24 hour mean minimum temperature (degC)
real(sp),    pointer, dimension(:)   :: tmax         ! 24 hour mean maximum temperature (degC)
real(sp),    pointer, dimension(:)   :: prec         ! 24 hour total precipitation (mm)
real(sp),    pointer, dimension(:)   :: cldf         ! 24 hour mean cloud cover (fraction)
real(sp),    pointer, dimension(:)   :: wind         ! 24 hour mean wind speed (m s-1)
real(sp),    pointer, dimension(:)   :: dayl         ! Daylength (h)
integer(i4), pointer, dimension(:)   :: sunrise      ! Current day sunrise hour (index of 24h hourly array)
integer(i4), pointer, dimension(:)   :: sunset       ! Current day sunset hour (index of 24h hourly array)
integer(i4), pointer, dimension(:)   :: sunrise_n    ! Next day sunrise hour (index of 24h hourly array)
integer(i4), pointer, dimension(:)   :: dayhour      ! Day time hours (h) --> from current day sunrise to sunset
integer(i4), pointer, dimension(:)   :: nighthour    ! Night time hours (h) --> from current day sunset to next day sunrise
real(sp),    pointer, dimension(:)   :: tmean        ! 24 hour mean temperature (degC)
real(sp),    pointer, dimension(:)   :: tday         ! Mean daytime temperature (degC)
real(sp),    pointer, dimension(:)   :: tnight       ! Mean nighttime temperature (degC)
real(sp),    pointer, dimension(:)   :: Ratm         ! Relative atmospheric pressure to sea-level (fraction)
real(sp),    pointer, dimension(:)   :: Ratm30       ! Relative atmospheric pressure to 30m above sea-level (fraction)
real(sp),    pointer, dimension(:)   :: Patm         ! Atmospheric pressure to sea-level (Pa)
real(sp),    pointer, dimension(:)   :: Patm30       ! Atmospheric pressure to 30m above sea-level (Pa)
real(sp),    pointer, dimension(:)   :: tdew         ! Dew point temperature (degC)
real(sp),    pointer, dimension(:)   :: Pjj          ! Precipitation equitability index for calculating PET
real(sp),    pointer, dimension(:)   :: rhum         ! Relative humidity (%)
real(sp),    pointer, dimension(:)   :: dsol         ! Solar declination angle (degree)
real(sp),    pointer, dimension(:)   :: srad         ! Downwelling surface shortwave radiation (kJ m-2 d-1)
real(sp),    pointer, dimension(:)   :: srad_dir     ! Direct beam downwelling shortwave raditaion (kJ m-2 d-1)
real(sp),    pointer, dimension(:)   :: srad_dif     ! Diffuse downwelling shortwave raditaion (kJ m-2 d-1)
real(sp),    pointer, dimension(:)   :: lrad         ! Upswelling surface longwave radiation (kJ m-2 d-1)
real(sp),    pointer, dimension(:)   :: dpet         ! Total potential evapotranspiration (mm)
real(sp),    pointer, dimension(:)   :: daet         ! Total actual evapotranspiration (mm)
real(sp),    pointer, dimension(:)   :: alpha        ! Ratio of daily AET/PET (fraction)
real(sp),    pointer, dimension(:)   :: vpd          ! Average daytime saturation vapor pressure deficit (Pa)
real(dp),    pointer, dimension(:,:) :: hprec        ! Hourly precipitation (mm)

logical,  pointer                 :: validcell
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
real(sp), pointer, dimension(:)   :: Tfield          ! Soil water volumetric content at field capacity (Psi = -33 kPa)   (fraction / m3 m-3)
real(sp), pointer, dimension(:)   :: Twilt           ! Soil water volumetric content at wilting point (Psi = -1500 kPa)   (fraction / m3 m-3)
real(sp), pointer, dimension(:)   :: Psat            ! Soil water matric potential at saturation (mm)
real(sp), pointer, dimension(:)   :: Bexp            ! Soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(sp), pointer, dimension(:)   :: Csolid          ! Soil solids volumetric heat capcity at layer midpoint (J m-3 K-1)
real(sp), pointer, dimension(:)   :: Ksolid          ! Soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
real(sp), pointer, dimension(:)   :: Kdry            ! Soil thermal conductivity of dry natural soil (W m-1 K-1)
real(sp), pointer, dimension(:)   :: Wliq            ! Soil liquid water content at layer midpoint (mm)
real(sp), pointer, dimension(:)   :: Wice            ! Soil ice content at layer midpoint (mm)
real(sp), pointer, dimension(:)   :: Tsoil           ! Soil temperature (K)
real(sp), pointer, dimension(:)   :: Tsoiln          ! Soil temperature for precious timestep (K)
real(sp), pointer, dimension(:)   :: Tliq            ! Soil volumetric water content (fraction / m3 m-3) / Fractional soil water content
real(sp), pointer, dimension(:)   :: Tice            ! Soil volumetric ice content (fraction / m3 m-3) / Fraction soil ice content
real(sp), pointer, dimension(:)   :: Psi             ! Soil water potential

real(sp), pointer                 :: toa_sw          ! Top of the atmosphere downwelling shortwave rad (kJ m-2 d-1)

real(sp), pointer, dimension(:,:) :: dwscal          ! Daily water stress factor (supply/demand ratio)
real(sp), pointer, dimension(:,:) :: dphen_t         ! Temperature based phenology of summergreen (proportion of leaf-on) (fraction)

logical,  pointer, dimension(:)   :: present         ! PFT present
logical,  pointer, dimension(:)   :: estab           ! PFT establishment
logical,  pointer, dimension(:)   :: survive         ! PFT survival
real(sp), pointer, dimension(:)   :: fpc_grid        ! Foilagrm  projective cover over grid (fraction)
real(sp), pointer, dimension(:)   :: fpc_ind         ! Foliage projective cover of individual (fraction)
real(sp), pointer, dimension(:)   :: fpc_inc         ! Foliage projective cover increment (fraction)
real(sp), pointer, dimension(:)   :: meanfpc
real(sp), pointer, dimension(:)   :: lm_ind          ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:)   :: rm_ind          ! Root carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:)   :: sm_ind          ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:)   :: hm_ind          ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:)   :: nind            ! PFT population
real(sp), pointer, dimension(:)   :: sla             ! Specific leaf area (m2 gC-1)
real(sp), pointer, dimension(:)   :: stemdiam        ! Tree stem diameter (m)
real(sp), pointer, dimension(:)   :: crownarea       ! Tree crownarea (m2)
real(sp), pointer, dimension(:)   :: height          ! Tree height (m)
real(sp), pointer, dimension(:)   :: lai_ind         ! Leaf area index of individual (m2 m-2)
real(sp), pointer, dimension(:)   :: litter_ag_fast  ! Fast above ground litter pool (gC m-2)
real(sp), pointer, dimension(:)   :: litter_ag_slow  ! Slow above ground litter pool (gC m-2)
real(sp), pointer, dimension(:)   :: litter_bg       ! Below ground litter pool (gC m-2)
logical,  pointer, dimension(:)   :: leafon          ! Leaf phenology on/off for photosynthesis
real(sp), pointer, dimension(:)   :: leafondays      ! Number of days since leaf phenology is on
real(sp), pointer, dimension(:)   :: leafoffdays     ! Number of days since leaf phenology is off

!-------------------------

real(sp), dimension(12) :: mtmp
real(sp), dimension(12) :: mpre

integer(i4) :: d

!-----------------------------------------------------------

lon       => sv(grid)%lon
lat       => sv(grid)%lat
georndst  => sv(grid)%georndst
validcell => sv(grid)%validcell

tmp => sv(grid)%monvars%tmp
dtr => sv(grid)%monvars%dtr
pre => sv(grid)%monvars%pre
wet => sv(grid)%monvars%wet
cld => sv(grid)%monvars%cld
wnd => sv(grid)%monvars%wnd
nd  => sv(grid)%monvars%nd

tmin      => sv(grid)%dayvars%tmin
tmax      => sv(grid)%dayvars%tmax
prec      => sv(grid)%dayvars%prec
cldf      => sv(grid)%dayvars%cldf
wind      => sv(grid)%dayvars%wind
dayl      => sv(grid)%dayvars%dayl
sunrise   => sv(grid)%dayvars%sunrise
sunset    => sv(grid)%dayvars%sunset
sunrise_n => sv(grid)%dayvars%sunrise_n
dayhour   => sv(grid)%dayvars%dayhour
nighthour => sv(grid)%dayvars%nighthour
tmean     => sv(grid)%dayvars%tmean
tday      => sv(grid)%dayvars%tday
tnight    => sv(grid)%dayvars%tnight
Ratm      => sv(grid)%dayvars%Ratm
Ratm30    => sv(grid)%dayvars%Ratm30
Patm      => sv(grid)%dayvars%Patm
Patm30    => sv(grid)%dayvars%Patm30
tdew      => sv(grid)%dayvars%tdew
Pjj       => sv(grid)%dayvars%Pjj
rhum      => sv(grid)%dayvars%rhum
dsol      => sv(grid)%dayvars%dsol
srad      => sv(grid)%dayvars%srad
srad_dir  => sv(grid)%dayvars%srad_dir
srad_dif  => sv(grid)%dayvars%srad_dif
lrad      => sv(grid)%dayvars%lrad
dpet      => sv(grid)%dayvars%dpet
daet      => sv(grid)%dayvars%daet
alpha     => sv(grid)%dayvars%alpha
vpd       => sv(grid)%dayvars%vpd
hprec     => sv(grid)%dayvars%hprec

dz        => sv(grid)%soilvars%dz
dzmm      => sv(grid)%soilvars%dzmm
zpos      => sv(grid)%soilvars%zpos
zposmm    => sv(grid)%soilvars%zposmm
zipos     => sv(grid)%soilvars%zipos
ziposmm   => sv(grid)%soilvars%ziposmm
sand      => sv(grid)%soilvars%sand
clay      => sv(grid)%soilvars%clay
cfvo      => sv(grid)%soilvars%cfvo
OrgM      => sv(grid)%soilvars%OrgM
rock      => sv(grid)%soilvars%rock
bulk      => sv(grid)%soilvars%bulk
Vsand     => sv(grid)%soilvars%Vsand
Vclay     => sv(grid)%soilvars%Vclay
Vsilt     => sv(grid)%soilvars%Vsilt
VOrgM     => sv(grid)%soilvars%VOrgM
awc       => sv(grid)%soilvars%awc
Ksat      => sv(grid)%soilvars%Ksat
Tsat      => sv(grid)%soilvars%Tsat
Tfield    => sv(grid)%soilvars%Tfield
Twilt     => sv(grid)%soilvars%Twilt
Psat      => sv(grid)%soilvars%Psat
Bexp      => sv(grid)%soilvars%Bexp
Csolid    => sv(grid)%soilvars%Csolid
Ksolid    => sv(grid)%soilvars%Ksolid
Kdry      => sv(grid)%soilvars%Kdry
Wliq      => sv(grid)%soilvars%Wliq
Wice      => sv(grid)%soilvars%Wice
Tsoil     => sv(grid)%soilvars%Tsoil
Tsoiln    => sv(grid)%soilvars%Tsoiln
Tliq      => sv(grid)%soilvars%Tliq
Tice      => sv(grid)%soilvars%Tice
Psi       => sv(grid)%soilvars%Psi

toa_sw    => sv(grid)%surfvars%toa_sw


dwscal  => sv(grid)%gppvars%dwscal
dphen_t => sv(grid)%gppvars%dphen_t

present        => sv(grid)%vegvars%present
estab          => sv(grid)%vegvars%estab
survive        => sv(grid)%vegvars%survive
fpc_grid       => sv(grid)%vegvars%fpc_grid
fpc_ind        => sv(grid)%vegvars%fpc_ind
fpc_inc        => sv(grid)%vegvars%fpc_inc
meanfpc        => sv(grid)%vegvars%meanfpc
lm_ind         => sv(grid)%vegvars%lm_ind
rm_ind         => sv(grid)%vegvars%rm_ind
sm_ind         => sv(grid)%vegvars%sm_ind
hm_ind         => sv(grid)%vegvars%hm_ind
nind           => sv(grid)%vegvars%nind
sla            => sv(grid)%vegvars%sla
stemdiam       => sv(grid)%vegvars%stemdiam
crownarea      => sv(grid)%vegvars%crownarea
height         => sv(grid)%vegvars%height
lai_ind        => sv(grid)%vegvars%lai_ind
litter_ag_fast => sv(grid)%vegvars%litter_ag_fast
litter_ag_slow => sv(grid)%vegvars%litter_ag_slow
litter_bg      => sv(grid)%vegvars%litter_bg
leafon         => sv(grid)%vegvars%leafon
leafondays     => sv(grid)%vegvars%leafondays
leafoffdays    => sv(grid)%vegvars%leafoffdays

elev => sv(grid)%topovars%elev

!-----------------------------------------------------------

if (yr == 1 .and. spinup) call initsoil(dz,dzmm,zpos,zposmm,zipos,ziposmm,sand,clay,cfvo,OrgM,rock,bulk, &
                           Vsand,Vclay,Vsilt,VOrgM,awc,Ksat,Tsat,Tfield,Twilt,Psat,Bexp, &
                           Csolid,Ksolid,Kdry,Wliq,Wice,Tsoil,Tsoiln,Tliq,Tice,Psi,validcell)

if (elev == -9999.) validcell = .false.

! if (yr == 1) call initsimplesoilLPJ(grid) ! Currently stored in waterbalancemod.f90

! Generate an initial grid-specific (geohash) random state for weathergen subroutine
! 'georndst' variable in randomdistmod
! Only need to be called at year 1 as the grid random state will be saved after gwgen() as pass onto the next year
if (yr == 1 .and. spinup) call initrndstate(lon,lat,georndst)

! Save current lon/lat of the grid in module variables 'clon' and 'clat' for later calculations
! NO LONGER IN USE (FOR NOW)
! call saveclonlat(grid)

! Copy 20 months (12 months +/- 4 months buffer) of monthly data for weathergen
call copymonvars(yr,grid,tmp,dtr,pre,wet,cld,wnd,nd)

! Generate daily met variables from monthly series (of original input variables)
call gwgen_new(georndst,tmp,dtr,pre,wet,cld,wnd,nd,tmin,tmax,prec,cldf,wind)

!-------------------------

mtmp = tmp(5:16)
mpre = pre(5:16)

!-------------------------

do d = 1, ndyear

  ! Subroutines to calculate other daily variables from tmin/tmax generated by gwgen()
  ! Calculate daylength, daytime and nighttime temperature
  call diurnaltemp(d,lat,tmin(d),tmin(d+1),tmax(d),tmean(d),tday(d),tnight(d),dayl(d),  &
                   sunrise(d),sunset(d),sunrise_n(d),dayhour(d),nighthour(d))

  ! Calculate atmospheric pressure variables based on elevation inputdata
  call elev_Ratm(elev,Ratm(d),Ratm30(d),Patm(d),Patm30(d),validcell)

  ! Calculate precipitation equitability index Pjj for tdew and pet in next subroutine
  call calcPjj(mtmp,mpre,Pjj(d))

  ! Calculate shortwave (direct + diffuse) radiation and longwave radiation (stable solution loop with pet)
  ! Calculate potential evapotranspiration and dewpoint temperature
  call radpet(d,lat,minval(mtmp),tmin(d),tmean(d),prec(d),cldf(d),Ratm(d),Pjj(d),tdew(d),dayl(d),  &
              dsol(d),toa_sw,srad(d),srad_dir(d),srad_dif(d),lrad(d),dpet(d))

  if (dayl(d) < 1.0)  dayl(d) = 1.0
  if (dayl(d) > 23.0) dayl(d) = 23.0

  ! Calculate estimated dewpoint from diurnal temperature range
  call calctdew(tmin(d),tmax(d),tmean(d),dayl(d),srad(d),dpet(d),sum(mpre),tdew(d))

  ! Calculate relative humidity
  call humidity(tmean(d),tdew(d),rhum(d))

  ! Calculate vapour pressure deficit
  call calcVPD(tmin(d),tday(d),vpd(d))

  ! Calculate actual evapotranspiration and AET/PET ration (alpha) using simple water bucket model (supply/demand from soil properties)
  ! ! call aet_alpha(yr,grid,d)

  ! Disaggregate 24 hour total prec into hourly series
  call hourlyprec(yr,grid,ndyear,d,prec(d),prec,tmean,rhum,hprec(d,:))

  ! if (tmin(d) > tmean(d)) print *, tdew(d),tmin(d),tmean(d),tmax(d),rhum(d),dayl(d),srad(d),dpet(d)

  ! Calculate fire danger meter
  ! call fireindex(rank,yr,grid,d)

  ! Print grid data if the process recieved the user-specified lon/lat grid
  ! if (lprint .AND. grid == gprint) call printgrid(info,grid,yr,d)

end do

!-------------------------

if (yr == 1 .and. spinup) call sapling(sla)
!
call bioclim(yr,grid,spinup,ndyear,mtmp,tmean(1:ndyear),estab,survive)

call summerphenology(ndyear,mtmp,tmean,dphen_t)

call establishment(yr,spinup,prec(1:ndyear),present,estab,survive,fpc_grid,fpc_ind,fpc_inc,lm_ind,rm_ind,sm_ind,hm_ind, &
                   nind,sla,stemdiam,height,crownarea,lai_ind,dwscal(d,:),leafon,leafondays,leafoffdays,         &
                   litter_ag_fast,litter_ag_slow,litter_bg)

end subroutine alvar_annual

!-------------------------------------------------------

end module alvarmod
