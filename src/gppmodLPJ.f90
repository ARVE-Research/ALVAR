module gppmod

! Module to calculate GPP coded for ALVAR by Leo Lai (Aug 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

subroutine gpp(year,grid,day)

! Based on combined methods in:
!     Haxeltine & Prentice (1996) BIOME3: An equilibrium terrestrial biosphere model based on ecophysiological constraints...  https://doi.org/10.1029/96GB02344
!     Wang et al. (2017) Towards a universal model for carbon dioxide uptake by plants, Nature Plnats. DOI:10.1038/s41477-017-0006-8
!
! Other references:
! 1. Photosynthetic photon flux density (PPFD)
!       Meek et al. (1984) A generalized relationship between PAR and SR. Agronomy J. https://doi.org/10.2134/agronj1984.00021962007600060018x
!       Jacovides et al. (2003) Global photosyntehtically active radiation and its relationship with global solar radition. Theor. Appl. Climatol. DOI 10.1007/s00704-002-0685-5
!       Yamashita & Yoshimura (2019) Estimation of global and diffuse PPFD. Remote Sens. http://dx.doi.org/10.3390/rs11080932
!       Mottus et al.(2011) Photosynthetically Active Radiation: Measurement and Modeling
! 2. Photorespiratory compensation point (gammastar)
!       Farquhar et al. (1980) A biochemical model of photosynthetic CO2 assimilation in leaves of C 3 species, Planta, 149, 78–90
!       Bernacchi et al. (2001) Improved temperature response functions for models of Rubisco-limited photosyn-thesis, Plant, Cell and Environment, 24, 253–259
! 3. Effective Michaelis–Menten coefficient of Rubisco
!       Bernacchi et al. (2001) Improved temperature response functions for models of Rubisco-limited photosyn-thesis, Plant, Cell and Environment, 24, 253–259

use utilitiesmod, only : getmonth
use pftparmod, only : npft,pftpar,tree,raingreen,c4
use statevarsmod, only : ndyear,dayvars,soilvars,gppvars,vegvars,topovars,gprint,lprint
use waterbalancemod, only : waterbalanceLPJ,w
use biome1mod,    only : biomevars
use soilstatemod, only : dzmm,dz

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

! Parameters
real(dp), parameter :: R        = 8.314      ! Universal gas constant (J mol-1 K-1)
real(dp), parameter :: phi0     = 1.02       ! Intrinsic quantum yield (g C mol-1) (Wang et al., 2017)
real(dp), parameter :: beta     = 240.       ! Constant in eq.3 (Wang et al. 2017)
real(dp), parameter :: beta_inv = 1./beta    ! Inverse of beta constant in eq.3 (Wang et al. 2017)
real(dp), parameter :: cstar    = 0.41       ! Unit carbon cost for maintenance of electron transport capacity (unitless) (Wang et al. 2017)
real(dp), parameter :: A        = -3.719     ! Constant in eq.13 for water viscosity (Wang et al. 2017)
real(dp), parameter :: B        = 580.       ! Constant in eq.13 for water viscosity (Wang et al. 2017)
real(dp), parameter :: C        = -138.      ! Constant in eq.13 for water viscosity (Wang et al. 2017)
real(dp), parameter :: sr2par   = 0.45       ! Conversion factor from solar irridiance (W m-2) to photosynthetically active radiation (PAR; 400-700 um) (W m-2) (Meek et al., 1984)
real(dp), parameter :: sr2ppfd  = 2.04       ! Conversion factor from solar irradiance (W m-2) to PPFD (umol m-2 s-1) (Meek et al., 1984)
real(dp), parameter :: coef_T   = 0.0545     ! Coefficient for temperature in eq.1 (Wang et al. 2017)
real(dp), parameter :: coef_D   = -0.5       ! Coefficient for vapor pressure deficit in eq.1 (Wang et al. 2017)
real(dp), parameter :: coef_z   = -0.0815    ! Coefficient for elevation in eq.1 (Wang et al. 2017)
real(dp), parameter :: intp_C   = 1.189      ! Intercept constant in eq.1 (Wang et al. 2017)
real(dp), parameter :: dhc      = 79430.     ! Activation energy for Kc (J mol-1) in eq.16 (Wang et al. 2017) --> from Bernacchi et al. 2001
real(dp), parameter :: dho      = 36380.     ! Activation energy for Ko (J mol-1) in eq.16 (Wang et al. 2017) --> from Bernacchi et al. 2001
real(dp), parameter :: kc25     = 39.97      ! k25 parameter for Kc variable (Bernacchi et al. 2001)
real(dp), parameter :: ko25     = 27480.     ! k25 parameter for Ko variable (Bernacchi et al. 2001)

!-------------------------
! Pointer variables
real(dp), pointer :: Ratm                 ! Relative atmospheric pressure to sea-level (fraction)
real(dp), pointer :: Patm                 ! Atmospheric pressure to sea-level (Pa)
real(dp), pointer :: vpd                  ! Average daytime saturation vapor pressure deficit (Pa)
real(dp), pointer :: dayl                 ! Daylength (h)
real(dp), pointer :: srad                 ! Downwelling surface shortwave radiation (kJ m-2 d-1)
real(dp), pointer :: srad_dir             ! Direct beam downwelling shortwave raditaion (kJ m-2 d-1)
real(dp), pointer :: srad_dif             ! Diffuse downwelling shortwave raditaion (kJ m-2 d-1)
real(dp), pointer :: tmean                ! 24 hour mean temperature (degC)
real(dp), pointer :: prec                 ! Daily precipitation (mm)
real(dp), pointer :: dpet                 ! Daily potential evapotranspiration (mm d-1)
real(dp), pointer :: daet
real(dp), pointer :: alpha

logical,  pointer, dimension(:) :: present              ! PFT present
real(dp), pointer, dimension(:) :: lambda               ! Actual ratio of lead internal to external CO2 (fraction)
real(dp), pointer, dimension(:) :: gammastar            ! Photorespiratory compensation point (gammastar in Wang et al. 2017) (Pa)
real(dp), pointer, dimension(:) :: gpp0                 ! Gross primary productivity under non-water stressed condition (g C m-2 d-1)
real(dp), pointer, dimension(:) :: gpp1                 ! Gross primary productivity under actual condition (g C m-2 d-1)
real(dp), pointer, dimension(:) :: gpp_tot              ! Total daily gross primary productivity (g C d-1)
real(dp), pointer, dimension(:) :: rd                   ! Daily leaf respiration (gC m-2 d-1) >> whole day include day + night
real(dp), pointer, dimension(:) :: dgp                  ! Optimal daily canopy conductance (mm m-2 s-1) ??? NEED CONFIRMATION ON UNIT in LPJ-LMFire
real(dp), pointer, dimension(:) :: dgc                  ! Actual daily canopy conductance (mm m-2 s-1) ??? NEED CONFIRMATION ON UNIT in LPJ-LMFire
real(dp), pointer, dimension(:) :: dwscal

real(dp), pointer, dimension(:) :: dphen
real(dp), pointer, dimension(:) :: dphen_t
real(dp), pointer, dimension(:) :: dphen_w

logical,  pointer, dimension(:) :: leafon
real(dp), pointer, dimension(:) :: leafondays
real(dp), pointer, dimension(:) :: leafoffdays

real(dp), pointer, dimension(:) :: fpc_grid

real(dp), pointer :: elev                 ! Elevation (m)
real(dp), pointer :: cellarea             ! Area of gridcell (m2)
real(dp), pointer :: areafrac             ! Ground area fraction in gridcell (fraction)

!-------------------------
! Local variables
real(dp) :: tmean_K               ! Daily mean temperature (K)

real(dp) :: xi                    ! Constant term in Wang et al. (2017) Eq. 27
real(dp) :: ci                    ! Optimal lead internal partial pressure of CO2 (Pa)
real(dp) :: ca                    ! Optimal external (ambient) partial pressure of CO2 (Pa)
real(dp) :: pi
real(dp) :: pa
real(dp) :: Kmm                   ! Effective Michaelis–Menten coefficient of Rubisco (Pa)
real(dp) :: Kc                    ! Michaelis constant for CO2 in eq.16 (Wang et al. 2017)
real(dp) :: Ko                    ! Michaelis constant for O2 in eq.16 (Wang et al. 2017)
real(dp) :: Po                    ! O2 partial pressure in eq.17 (Pa) (Wang et al. 2017)
real(dp) :: visco                 ! Viscosity of water relative to 25C / 298.15K (fraction)
real(dp) :: par                   ! Photosynthetically active radiation (W m-2)
real(dp) :: fpar                  ! Fraciton of PAR intercepted by foliage (fraction)
real(dp) :: apar                  ! Absorbed photosynthetically active radiation (W m-2)
real(dp) :: ppfd                  ! Photosynthetic photo flux density (mol m-2 s-1)
real(dp) :: fppfd                 ! Fraciton of PPFD intercepted by foliage (fraction)
real(dp) :: appfd                 ! Absorbed (by leaf) photosynthetic photo flux density (mol m-2 s-1)
real(dp) :: lambdam

real(dp) :: rootprop              ! Root proportion in soil layer (fraction)
real(dp) :: gminp                 ! Minimum canopy conductance rate (mm m-2 s-1)

real(dp) :: tau
real(dp) :: tsecs
real(dp) :: supply                ! Supply of soil moisture (mm d-1)
real(dp) :: demand                ! Demand for soil moisture (mm d-1)
real(dp) :: vm0                   ! Optimal rubisco activity (gC m-2 d-1)
real(dp) :: rd0                   ! Optimal leaf respiration rate (gC m-2 d-1)
real(dp) :: and0                  ! Optimal daily net photosynthesis (gC m-2 d-1)
real(dp) :: adt0                  ! Optimal total daytime net photosynthesis (gC m-2 d-1)
real(dp) :: adtmm0                ! Optimal total daytime net photosynthesis in mm (mm m-2 d-1)
real(dp) :: vm1                   ! Actual leaf rubisco activity (gC m-2 d-1)
real(dp) :: adtmm1                ! Actual total daytime net photosynthesis in mm (mm m-2 d-1)
real(dp) :: adtmm2
real(dp) :: je                    ! APAR-limited photosynthesis rate (molC m-2 h-1)
real(dp) :: jc                    ! Rubisco-activity-limited photosynthesis rate (molC m-2 h-1)
real(dp) :: agd0
real(dp) :: agd1

!-------------------------
! Temperature inhibition variables
real(dp) :: inhibx1
real(dp) :: inhibx2
real(dp) :: inhibx3
real(dp) :: inhibx4
real(dp) :: k1
real(dp) :: k2
real(dp) :: k3
real(dp) :: low
real(dp) :: high
real(dp) :: tstress

! Rubisco capcity variables
real(dp) :: c1
real(dp) :: c2
real(dp) :: b0
real(dp) :: t0
real(dp) :: s
real(dp) :: sigma

! Bisection root finding variables
real(dp) :: gpd
real(dp) :: epsilon
real(dp) :: x1
real(dp) :: x2
real(dp) :: rtbis
real(dp) :: dx
real(dp) :: fmid
real(dp) :: xmid

real(dp) :: longevity
real(dp) :: minwscal
real(dp) :: dwscal_ytd
real(dp) :: drunoff_surf
real(dp) :: drunoff_drain

integer :: pft
integer :: it

!-------------------------

Ratm     => dayvars(grid,day)%Ratm
Patm     => dayvars(grid,day)%Patm
vpd      => dayvars(grid,day)%vpd
dayl     => dayvars(grid,day)%dayl
srad     => dayvars(grid,day)%srad
srad_dir => dayvars(grid,day)%srad_dir
srad_dif => dayvars(grid,day)%srad_dif
tmean    => dayvars(grid,day)%tmean
prec     => dayvars(grid,day)%prec
dpet     => dayvars(grid,day)%dpet
daet     => dayvars(grid,day)%daet
alpha    => dayvars(grid,day)%alpha

elev     => topovars(grid)%elev
cellarea => topovars(grid)%cellarea
areafrac => topovars(grid)%areafrac

gammastar=> gppvars(grid,day)%gammastar
lambda   => gppvars(grid,day)%lambda
gpp0     => gppvars(grid,day)%gpp0
gpp1     => gppvars(grid,day)%gpp
gpp_tot  => gppvars(grid,day)%gpp_tot
rd       => gppvars(grid,day)%rd
dgp      => gppvars(grid,day)%dgp
dgc      => gppvars(grid,day)%dgc
dwscal   => gppvars(grid,day)%dwscal

present  => vegvars(grid)%present
dphen    => vegvars(grid)%dphen(day,:)
dphen_t  => vegvars(grid)%dphen_t(day,:)
dphen_w  => vegvars(grid)%dphen_w(day,:)

leafon      => vegvars(grid)%leafon
leafondays  => vegvars(grid)%leafondays
leafoffdays => vegvars(grid)%leafoffdays

fpc_grid => vegvars(grid)%fpc_grid

!-------------------------

if (areafrac <= 0.) then
  gpp1 = -9999.
  gpp_tot = 0.
  return
end if

tmean_K = tmean + Tfreeze

if (elev < 0.) elev = 0.
if (areafrac < 0.) areafrac = 0.

!-------------------------

ca = 4.00e-4
! ca = 2.88e-4 + 0.009e-4 * year
par = srad * 1000. * 0.5     ! Convert from kJ m-2 d-1 to J m-2 d-1 and divide by half to get PAR
tsecs = 3600. * dayl

!---------------------------------------------------------------------
! Calculate photosynthesis (agd0) and canopy conductance (dgp) under non-water-stressed conditions

do pft = 1, npft

  if (present(pft)) then

    !------

    inhibx1 = pftpar(22,pft)
    inhibx2 = pftpar(23,pft)
    inhibx3 = pftpar(24,pft)
    inhibx4 = pftpar(25,pft)

    gminp   = pftpar(4,pft) * fpc_grid(pft)
    lambdam = pftpar(26,pft)

    !------

    fpar = fpc_grid(pft)

    call photosynthesis(pft,ca,dayl,fpar,lambdam,c4(pft),par,tmean,inhibx1,inhibx2,inhibx3,inhibx4,adtmm0,agd0,rd(pft))

    if (tsecs > 0.) then
      dgp(pft) = (((1.6 * adtmm0) / (ca * (1. - lambdam))) / tsecs) + gminp
    else
      dgp(pft) = 0.
    end if

    ! Save optimal gpp value to statevars
    gpp0(pft) = agd0

  end if

end do

!---------------------------------------------------------------------
! Calculate water-stressed leaf phenology based on water scalr from previous day
! then calculate water balance for the day
! Also obtain actual canopy conductance (dgc) from water balance subroutine
do pft = 1, npft

  if (present(pft)) then

    minwscal  = pftpar(3,pft)
    longevity = pftpar(7,pft)

    if (day == 1) then
      dwscal_ytd = gppvars(grid,365)%dwscal(pft)
    else
      dwscal_ytd = gppvars(grid,day-1)%dwscal(pft)
    end if

    !------

    if (dwscal_ytd > minwscal .and. leafon(pft)) then

      dphen_w(pft) = 1.
      dphen(pft) = dphen_t(pft)
      leafondays(pft) = leafondays(pft) + 1

    else

      dphen_w(pft) = 0.
      dphen(pft) = 0.

    end if

    !------
    ! Stop deciduous vegetation (raingreen) behaving like evergreen when climate permits

    if (raingreen(pft) .and. tree(pft)) then

        if (real(leafondays(pft)) >= (365. * longevity)) then

          leafon(pft)      = .false.
          leafoffdays(pft) = leafoffdays(pft) + 1

          if (real(leafoffdays(pft)) >= (365. * longevity)) then

            leafoffdays(pft) = 0
            leafondays(pft)  = 0
            leafon(pft)      = .true.

          end if

        end if

      end if

      !------

  end if

end do

!------
! Update current day soilwater balance based on precipitation and soil properties
! NOTE: Currently in LPJ-LMFire method, where all related soil variables are self-contained in waterbalancemod.f90
!       Preferably switch to the 6 layer soil model coded in hydrologymod.f90 some time soon
call waterbalanceLPJ(grid,day,present,dgp,dpet,dphen,dgc,prec,drunoff_drain,drunoff_surf,dwscal,daet,fpc_grid)


!---------------------------------------------------------------------
! Calculate photosynthesis in water stressed condition by considering actual canopy conductance (gpd)
do pft = 1, npft

  if (present(pft)) then

    !------

    inhibx1 = pftpar(22,pft)
    inhibx2 = pftpar(23,pft)
    inhibx3 = pftpar(24,pft)
    inhibx4 = pftpar(25,pft)

    gminp   = pftpar(4,pft) * fpc_grid(pft)
    lambdam = pftpar(26,pft)

    !------

    fpar = fpc_grid(pft)

    !------

    gpd = tsecs * (dgc(pft) - gminp)      ! Convert canopy conductance from mm m-2 s-1 to mm m-2 d-1

    !-------------------------
    ! Bisection root finding method to solve Eq. 2 and Eq. 18 simultaneuously (Haxeltine & Prentice, 1996)
    ! To find actual intercellular to ambient partial pressure of CO2 (chi, or lamba in LPJ)
    epsilon = 0.05

    if (gpd > 1.e-5) then     ! Significant canopy conductance

      ! Initiate bracket for numerical solution

      x1 = 0.02               ! Minimum bracket of the root
      x2 = lambdam + 0.05      ! Maximum bracket of the root = optimal ratio
      rtbis = x1              ! Root of the bisection
      dx = x2 - x1

      it = 0                  ! Number of tries towards solution

      fmid = epsilon + 1.

      do  ! Bisection root finding

        it  = it + 1
        dx  = dx * 0.5
        xmid = rtbis + dx

        ! Calculate total daytime photosynthesis implied by canopy conductance from dgc calculation and
        ! current guess for lambda (xmid).
        ! Units are mm m-2 d-1 (mm come from gpd value, mm d-1)
        ! Reverse of Eq. 18 (Haxeltine & Prentice, 1996)
        adtmm1 = gpd / 1.6 * ca * (1. - xmid)

        ! Evaluate fmid at the point lambda = xmid fmid will be an increasing function with xmid,
        ! with a solution (fmid=0) between x1 and x2
        call photosynthesis(pft,ca,dayl,fpar,xmid,c4(pft),par,tmean,inhibx1,inhibx2,inhibx3,inhibx4,adtmm2,agd1,rd(pft))
        ! call photosynthesis(pft,ca,dayl,fpar,lambdam,c4(pft),par,tmean,inhibx1,inhibx2,inhibx3,inhibx4,adtmm2,agd1,rd(pft))

        fmid = adtmm2 - adtmm1

        ! if (lprint .and. grid==gprint .and. pft==1) write(0,*) fmid, adtmm2, adtmm1, agd1, gpp0(pft)

        if (fmid < 0.) rtbis = xmid

        ! Exit if solution found or > 10 iterations needed without solution
        if (abs(fmid) < epsilon) exit

        if (it > 10) then
          if (lprint .and. grid==gprint .and. pft==1) write(0,*) 'No photosynthesis solution for pft: ', pft, it
          exit
        end if

      end do ! End of bisection

    else  ! Infinitesimal canopy conductance

      rd(pft)   = 0.
      agd1  = 0.
      xmid = 0.

    end if  ! Canopy conductance IF condition

    !---------------------------------------------------------------------
    ! Save values into state variables
    ! gpp0(pft)   = agd0
    gpp1(pft)   = agd1
    lambda(pft) = xmid
    alpha       = min(daet/dpet, 1.)

    gpp_tot(pft) = gpp1(pft) * cellarea * areafrac

  else

    gpp0(pft)   = 0.
    gpp1(pft)   = 0.
    lambda(pft) = 0.
    alpha       = min(daet/dpet, 1.)

    gpp_tot(pft) = gpp1(pft) * cellarea * areafrac

  end if ! Present IF condition

end do ! PFT loop


end subroutine gpp

!---------------------------------------------------------------------

subroutine photosynthesis(pft,ca,dayl,fpar,lambda,c4,par,temp,x1,x2,x3,x4,adtmm,agd,rd)

use parametersmod, only : sp,dp
use pftparmod,     only : npft,pftpar

implicit none

! Parameters copied from LPJ-LMFire (photosynthesismod)
real(dp), parameter :: alphaa    =    0.5     ! fraction of PAR assimilated at ecosystem level
real(dp), parameter :: alphac3   =    0.08    ! intrinsic quantum efficiency of CO2 uptake in C3 plants
real(dp), parameter :: alphac4   =    0.053   ! C4 intrinsic quantum efficiency
real(dp), parameter :: bc3       =    0.015   ! leaf respiration as fraction of Vmax for C3 plants
real(dp), parameter :: bc4       =    0.02    ! leaf respiration as fraction of vmax for C4 plants
real(dp), parameter :: cmass     =   12.      ! atomic mass of carbon
real(dp), parameter :: cq        =    4.6e-6  ! conversion factor for solar radiation at 550 nm from J/m2 to E/m2 (E=mol quanta)
real(dp), parameter :: e0        =  308.56    ! parameter in Arrhenius temp response function
real(dp), parameter :: kc25      =   30.      ! value of kc at 25 deg C
real(dp), parameter :: ko25      =    3.e4    ! value of ko at 25 deg C
real(dp), parameter :: lambdamc3 =    0.8     ! optimal (maximum?) ci:ca ratio for C3 plants
real(dp), parameter :: lambdamc4 =    0.4     ! optimal ci:ca ratio for c4 plantsfpc_gr
real(dp), parameter :: m         =   25.      ! corresponds to parameter p in Eqn 28, Haxeltine & Prentice 1996
real(dp), parameter :: n0        =    7.15    ! leaf N concentration (mg/g) not involved in photosynthesis
real(dp), parameter :: p         =    1.e5    ! atmospheric pressure (Pa)
real(dp), parameter :: po2       =   20.9e3   ! O2 partial pressure (Pa)
real(dp), parameter :: q10kc     =    2.1     ! q10 for temperature-sensitive parameter kc
real(dp), parameter :: q10ko     =    1.2     ! q10 for temperature-sensitive parameter ko
real(dp), parameter :: q10tau    =    0.57    ! q10 for temperature-sensitive parameter tau
real(dp), parameter :: t0c3      =  250.      ! base temperature (K) in Arrhenius temperature response function for C3 plants
real(dp), parameter :: t0c4      =  260.      ! base temperature in Arrhenius func for C4 plants
real(dp), parameter :: tau25     = 2600.      ! value of tau at 25 deg C
real(dp), parameter :: theta     =    0.7     ! colimitation (shape) parameter
real(dp), parameter :: tk25      =  298.15    ! 25 deg C in Kelvin
real(dp), parameter :: tmc3      =   45.      ! maximum temperature for C3 photosynthesis
real(dp), parameter :: tmc4      =   55.      ! maximum temperature for C4 photosynthesis

! Subroutine arguments
! logical,  intent(in)  :: c4
! integer,  intent(in)  :: pfti
integer,  intent(in)  :: pft
real(dp), intent(in)  :: ca
real(dp), intent(in)  :: dayl
real(dp), intent(in)  :: fpar
real(dp), intent(in)  :: lambda
logical,  intent(in)  :: c4
real(dp), intent(in)  :: par
real(dp), intent(in)  :: temp
real(dp), intent(in)  :: x1
real(dp), intent(in)  :: x2
real(dp), intent(in)  :: x3
real(dp), intent(in)  :: x4
real(dp), intent(out) :: adtmm
real(dp), intent(out) :: agd
real(dp), intent(out) :: rd

! Local variables
real(dp) :: adt
real(dp) :: and
real(dp) :: apar
real(dp) :: b
real(dp) :: c1
real(dp) :: c2
real(dp) :: gammastar
real(dp) :: high
real(dp) :: jc
real(dp) :: je
real(dp) :: k1
real(dp) :: k2
real(dp) :: k3
real(dp) :: kc
real(dp) :: ko
real(dp) :: lambdam
real(dp) :: low
real(dp) :: pa
real(dp) :: phipi
real(dp) :: pi
real(dp) :: s
real(dp) :: sigma
real(dp) :: t0
real(dp) :: tau
real(dp) :: tstress
real(dp) :: vm

!-------------------------
! Return withtou perfomring calculations if daylength is too short

if (dayl < 0.01) then

  agd   = 0.
  adtmm = 0.
  rd    = 0.

  return

end if

!-------------------------
! Calculate APAR in J/m2/day (Eq. 4 in Hazeltine & Prentice, 1996)
! alphaa = scaling factor for absorbed PAR at ecosystem

apar = par * fpar * alphaa

lambdam = pftpar(26,pft) ! Get optimal Ci/Ca ratio for PFT

!-------------------------
! Calculate temperature inhibition function on photosynthesis
! High-temperature inhibition modelled conservatively as a step function prohibiting photosynthesis above 55 deg C (Table 3.7, Larcher 1983)
! NOTE: C3 plant assupmtion for now

if (temp < x4) then

   k1  = 2. * alog((1. / 0.99) - 1.) / (x1 - x2)
   k2  = (x1 + x2) / 2.
   low = 1. / (1. + exp(k1 * (k2 - temp)))

   k3   = alog(0.99 / 0.01) / (x4 - x3)
   high = 1. - 0.01 * exp(k3 * (temp - x3))

   tstress = low * high

else    ! High temperature stress on photosynthesis, tstress = 0.

   tstress = 0.

end if

if (tstress < 1.e-2) tstress = 0.

!-------------------------
! First calculate catalytic capacity of rubisco, Vm, assuming optimal (non-water-stressed) value for lambda, i.e. lambdamc3

if (c4) then   ! C4 photosynthesis (grass)

  c1 = tstress * alphac4

  ! High-temperature inhibition modelled conservatively as a step function
  ! prohibiting photosynthesis above 55 deg C (Table 3.7, Larcher 1983)

  if (temp > tmc4) c1 = 0.

  c2 = 1.

  b  = bc4    !Choose C4 value of b for Eqn 10, Haxeltine & Prentice 1996
  t0 = t0c4   !base temperature for temperature response of rubisco

else ! C3 photosynthesis

  !-------------------------
  ! Recalculate parameters from LPJ paper method (Leo Lai Oct 2021)
  ko  = ko25  * q10ko **((temp - 25.) / 10.)  !Michaelis constant of rubisco for O2
  kc  = kc25  * q10kc **((temp - 25.) / 10.)  !Michaelis constant for CO2
  tau = tau25 * q10tau**((temp - 25.) / 10.)  !CO2/O2 specificity ratio

  !-------------------------
  ! Calculate CO2 compensation point gammastar (CO2 partial pressure; Pa)
  gammastar = po2 / (2. * tau)                                            ! Eq. 8 (Haxeltine & Prentice, 1996)

  !-------------------------
  ! Non-water-stressed intercellular CO2 partiual pressure (Pa)
  pa = ca * p

  pi = lambdam * pa

  !-------------------------
  ! Calculation of optimal rubisco capacity, Vm (gC/m2/day) from BIOME3 model (Haxeltine & Prentice, 1996)
  c1 = tstress * alphac3 * ((pi - gammastar) / (pi + 2. * gammastar))     ! Eq. 4 (Haxeltine & Prentice, 1996)

    if (temp > tmc3) c1 = 0.

  c2 = (pi - gammastar) / (pi + kc * (1.0 + po2 / ko))                     ! Eq. 6 (Haxeltine & Prentice, 1996)

  b  = bc3
  t0 = t0c3

end if

! Calculate catalytic capacity of rubisco, Vm, assuming optimal (non-water-stressed) value for lambda, i.e. lambdamc3
s = (24. / dayl) * b                                                    ! Eq. 13 (Haxeltine & Prentice, 1996)

sigma = sqrt(max(0., 1. - (c2 - s) / (c2 - theta * s)))                 ! Eq. 12 (Haxeltine & Prentice, 1996)

vm  = (1. / b) * (c1 / c2) * ((2. * theta - 1.) * s - &                 ! Eq. 11 (Haxeltine & Prentice, 1996)
     (2. * theta * s - c2) * sigma) * apar * cmass * cq

!-------------------------
! Use the Vm value above to calculate actual photosynthesis and actual lambda for actual pi CO2 partial pressure (Pa)
! Intercellular CO2 partial pressure in Pa

if (.not.c4) then

  pi = lambda * pa                                                        ! Eq. 7 (Haxeltine & Prentice, 1996)

  ! Recalculation of C1C3, C2C3 with actual pi

  c1 = tstress * alphac3 * ((pi - gammastar) / (pi + 2. * gammastar))

  if (temp > tmc3) c1 = 0.  !High-temperature inhibition

  c2 = (pi - gammastar) / (pi + kc * (1. + po2 / ko))

else

  ! Parameter accounting for effect of reduced intercellular CO2 concentration on photosynthesis, Phipi
  ! Eqn 14,16, Haxeltine & Prentice 1996; Fig 1b, Collatz et al 1992

  phipi = min(lambda / lambdam, 1.)

  c1 = tstress * phipi * alphac4

  if (temp > tmc4) c1 = 0.  ! High-temperature inhibition

end if


!-------------------------
! Calculation of APAR-limited photosynthesis rate Je, molC m-2 -h;              ! Eq. 3 (Haxeltine & Prentice, 1996)
je = c1 * apar * cmass * cq / dayl

! Calculation of rubisco-activity-limited photosynthesis rate Jc, molC m-2 -h   ! Eq. 5 (Haxeltine & Prentice, 1996)
jc = c2 * vm / 24.

if (je < 1.e-10 .or. jc <= 1.e-10) then

  agd = 0.0

else

  ! Calculation of daily gross photosynthesis, Agd, gC/m2/day; Eq. 2 (Haxeltine & Prentice, 1996)
  ! NOTE: there is an error in this equation in the above paper (missing theta in 4*theta*je*jc term) which is fixed here (LPJ-LMFire comment)

  agd = (je + jc - sqrt((je + jc)**2. - 4. * theta * je * jc)) / (2. * theta) * dayl

end if

! Optimal daily leaf respiration, Rd, gC m-2 d-1                                ! Eq. 10 (Haxeltine & Prentice, 1996)
rd = b * vm

! Optimal daily net photosynthesis (at leaf level),And, gC m-2 d-1 >> using the optimal condition estimated gpp0 from (Wang et al., 2017)
and = agd - rd

! Optimal total daytime net photosynthesis, Adt, gC m-2 d-1                     ! Eq. 19 (Haxeltine & Prentice, 1996)
adt = and + (1. - dayl / 24.) * rd

! Convert adt from gC m-2 d-1 to mm m-2 d-1 using ideal gas equation
adtmm = ((adt / cmass) * 8.314 * (temp + 273.3) / p) * 1000.


end subroutine photosynthesis

!---------------------------------------------------------------------
!
! subroutine leafcarbon(grid,day)
!
! use statevarsmod, only : vegvars
!
! implicit none
!
! integer(i4), intent(in) :: grid
! integer(i4), intent(in) :: day
!
! real(dp), parameter, dimension(17) :: sla0 = [0.012    &       ! 1 = Tropical evergreen / Broadlead evergreen tree (tropical)
!                                               ,0.012    &       ! 2 = Tropical raingreen / Broadleaf evergreen tree (tropical)
!                                               ,0.030    &       ! 3 = Tropical deciduous / Broadleaf deciduous tree (tropical)
!                                               ,0.012    &       ! 4 = Temperate evergreen (warm mixed) / Broadleaf evergreen tree (temperate)
!                                               ,0.030    &       ! 5 = Temperate deciduous / Broadleaf deciduous tree (temperate)
!                                               ,0.030    &       ! 6 = Cool mixed / Broadleaf deciduous tree (temperate)
!                                               ,0.010    &       ! 7 = Cool conifer / Needleleaf evergreen (temperate)
!                                               ,0.008    &       ! 8 = Cold evergreen / Needleleaf evergreen tree (boreal)
!                                               ,0.008    &       ! 9 = Cold mixed / Needleleaf evergreen tree (boreal)
!                                               ,0.024    &       ! 10 = Cold deciduous / Needleleaf deciduous tree
!                                               ,0.050    &       ! 11 = Xerophytic / C4 grass
!                                               ,0.050    &       ! 12 = Warm grass / C4 grass
!                                               ,0.050    &       ! 13 = Cool grass / C3 nonarctic grass
!                                               ,0.000    &       ! 14 = Tundra
!                                               ,0.000    &       ! 15 = Hot desert
!                                               ,0.000    &       ! 16 = Cold desert
!                                               ,0.000    ]       ! 17 = Barren
!
! real(dp), parameter, dimension(17) :: m    = [0.0015    &       ! 1 = Tropical evergreen / Broadlead evergreen tree (tropical)
!                                               ,0.0015    &       ! 2 = Tropical raingreen / Broadleaf evergreen tree (tropical)
!                                               ,0.0040    &       ! 3 = Tropical deciduous / Broadleaf deciduous tree (tropical)
!                                               ,0.0015    &       ! 4 = Temperate evergreen (warm mixed) / Broadleaf evergreen tree (temperate)
!                                               ,0.0040    &       ! 5 = Temperate deciduous / Broadleaf deciduous tree (temperate)
!                                               ,0.0040    &       ! 6 = Cool mixed / Broadleaf deciduous tree (temperate)
!                                               ,0.0013    &       ! 7 = Cool conifer / Needleleaf evergreen (temperate)
!                                               ,0.0010    &       ! 8 = Cold evergreen / Needleleaf evergreen tree (boreal)
!                                               ,0.0010    &       ! 9 = Cold mixed / Needleleaf evergreen tree (boreal)
!                                               ,0.0030    &       ! 10 = Cold deciduous / Needleleaf deciduous tree
!                                               ,0.0000    &       ! 11 = Xerophytic / C4 grass
!                                               ,0.0000    &       ! 12 = Warm grass / C4 grass
!                                               ,0.0000    &       ! 13 = Cool grass / C3 nonarctic grass
!                                               ,0.0000    &       ! 14 = Tundra
!                                               ,0.0000    &       ! 15 = Hot desert
!                                               ,0.0000    &       ! 16 = Cold desert
!                                               ,0.0000    ]       ! 17 = Barren
!
! integer(i4), pointer :: biome
! real(dp), pointer :: C_leaf
! real(dp), pointer :: sla
! real(dp), pointer :: lai
!
! !-------------------------
!
! C_leaf => vegvars(grid,day)%C_leaf
! sla => vegvars(grid,day)%sla
! lai => vegvars(grid,day)%lai
! biome => vegvars(grid,day)%biome
!
! !-------------------------
!
! if (biome /= missing_i2 .and. biome < 11 .and. lai /= -9999.) then
!
!   C_leaf = log(m(biome) * lai + sla0(biome) - log(sla0(biome))) / m(biome)
!
! end if
!
!
!
! end subroutine leafcarbon
!
! !---------------------------------------------------------------------
!
! subroutine lai(grid,day)
!
! use statevarsmod, only : vegvars,gprint,lprint
!
! implicit none
!
! integer(i4), intent(in) :: grid
! integer(i4), intent(in) :: day
!
! real(dp), parameter, dimension(17) :: sla0 = [0.012    &       ! 1 = Tropical evergreen / Broadlead evergreen tree (tropical)
!                                               ,0.012    &       ! 2 = Tropical raingreen / Broadleaf evergreen tree (tropical)
!                                               ,0.030    &       ! 3 = Tropical deciduous / Broadleaf deciduous tree (tropical)
!                                               ,0.012    &       ! 4 = Temperate evergreen (warm mixed) / Broadleaf evergreen tree (temperate)
!                                               ,0.030    &       ! 5 = Temperate deciduous / Broadleaf deciduous tree (temperate)
!                                               ,0.030    &       ! 6 = Cool mixed / Broadleaf deciduous tree (temperate)
!                                               ,0.010    &       ! 7 = Cool conifer / Needleleaf evergreen (temperate)
!                                               ,0.008    &       ! 8 = Cold evergreen / Needleleaf evergreen tree (boreal)
!                                               ,0.008    &       ! 9 = Cold mixed / Needleleaf evergreen tree (boreal)
!                                               ,0.024    &       ! 10 = Cold deciduous / Needleleaf deciduous tree
!                                               ,0.050    &       ! 11 = Xerophytic / C4 grass
!                                               ,0.050    &       ! 12 = Warm grass / C4 grass
!                                               ,0.050    &       ! 13 = Cool grass / C3 nonarctic grass
!                                               ,0.000    &       ! 14 = Tundra
!                                               ,0.000    &       ! 15 = Hot desert
!                                               ,0.000    &       ! 16 = Cold desert
!                                               ,0.000    ]       ! 17 = Barren
!
! real(dp), parameter, dimension(17) :: m    = [0.0015    &       ! 1 = Tropical evergreen / Broadlead evergreen tree (tropical)
!                                               ,0.0015    &       ! 2 = Tropical raingreen / Broadleaf evergreen tree (tropical)
!                                               ,0.0040    &       ! 3 = Tropical deciduous / Broadleaf deciduous tree (tropical)
!                                               ,0.0015    &       ! 4 = Temperate evergreen (warm mixed) / Broadleaf evergreen tree (temperate)
!                                               ,0.0040    &       ! 5 = Temperate deciduous / Broadleaf deciduous tree (temperate)
!                                               ,0.0040    &       ! 6 = Cool mixed / Broadleaf deciduous tree (temperate)
!                                               ,0.0013    &       ! 7 = Cool conifer / Needleleaf evergreen (temperate)
!                                               ,0.0010    &       ! 8 = Cold evergreen / Needleleaf evergreen tree (boreal)
!                                               ,0.0010    &       ! 9 = Cold mixed / Needleleaf evergreen tree (boreal)
!                                               ,0.0030    &       ! 10 = Cold deciduous / Needleleaf deciduous tree
!                                               ,0.0000    &       ! 11 = Xerophytic / C4 grass
!                                               ,0.0000    &       ! 12 = Warm grass / C4 grass
!                                               ,0.0000    &       ! 13 = Cool grass / C3 nonarctic grass
!                                               ,0.0000    &       ! 14 = Tundra
!                                               ,0.0000    &       ! 15 = Hot desert
!                                               ,0.0000    &       ! 16 = Cold desert
!                                               ,0.0000    ]       ! 17 = Barren
!
! integer(i4), pointer :: biome
! real(dp), pointer :: C_leaf
! real(dp), pointer :: sla
! real(dp), pointer :: lai0
!
! real(dp) :: lai_sun
! real(dp) :: lai_sha
!
! real(dp), parameter :: kb = 0.5       ! Extinction factors, should be function of zenith angle and leaf angle, but constant FOR NOW (Dai et al., 2004)
!
!
! !-------------------------
!
! C_leaf => vegvars(grid,day)%C_leaf
! sla => vegvars(grid,day)%sla
! lai0 => vegvars(grid,day)%lai
! biome => vegvars(grid,day)%biome
!
! !-------------------------
!
! if (biome /= missing_i2 .and. biome < 11 .and. lai0 /= -9999.) then
!
!   lai0 = sla0(biome) * exp(m(biome) * C_leaf - 1.0) / m(biome)
!
! end if
!
! lai_sun = (1 - exp(-kb * lai0)) / kb
! lai_sha = lai0 - lai_sun
!
! if(grid == gprint .and. lprint) print *, day, lai0, lai_sun,lai_sha, C_leaf
!
! end subroutine lai
!
! !---------------------------------------------------------------------
!
! subroutine gpp_new(year,grid,day)
!
! ! Based on combined methods in:
! !     Haxeltine & Prentice (1996) BIOME3: An equilibrium terrestrial biosphere model based on ecophysiological constraints...  https://doi.org/10.1029/96GB02344
! !     Wang et al. (2017) Towards a universal model for carbon dioxide uptake by plants, Nature Plnats. DOI:10.1038/s41477-017-0006-8
! !
! ! Other references:
! ! 1. Photosynthetic photon flux density (PPFD)
! !       Meek et al. (1984) A generalized relationship between PAR and SR. Agronomy J. https://doi.org/10.2134/agronj1984.00021962007600060018x
! !       Jacovides et al. (2003) Global photosyntehtically active radiation and its relationship with global solar radition. Theor. Appl. Climatol. DOI 10.1007/s00704-002-0685-5
! !       Yamashita & Yoshimura (2019) Estimation of global and diffuse PPFD. Remote Sens. http://dx.doi.org/10.3390/rs11080932
! !       Mottus et al.(2011) Photosynthetically Active Radiation: Measurement and Modeling
! ! 2. Photorespiratory compensation point (gammastar)
! !       Farquhar et al. (1980) A biochemical model of photosynthetic CO2 assimilation in leaves of C 3 species, Planta, 149, 78–90
! !       Bernacchi et al. (2001) Improved temperature response functions for models of Rubisco-limited photosyn-thesis, Plant, Cell and Environment, 24, 253–259
! ! 3. Effective Michaelis–Menten coefficient of Rubisco
! !       Bernacchi et al. (2001) Improved temperature response functions for models of Rubisco-limited photosyn-thesis, Plant, Cell and Environment, 24, 253–259
!
! use utilitiesmod, only : getmonth
! use statevarsmod,   only : ndyear,dayvars,soilvars,vegvars,topovars,gprint,lprint
! use biome1mod,    only : biomevars
! use soilstatemod, only : dzmm,dz
!
! integer(i4), intent(in) :: year
! integer(i4), intent(in) :: grid
! integer(i4), intent(in) :: day
!
! ! Parameters
! real(dp), parameter :: R        = 8.314      ! Universal gas constant (J mol-1 K-1)
! real(dp), parameter :: phi0     = 1.02       ! Intrinsic quantum yield (g C mol-1) (Wang et al., 2017)
! real(dp), parameter :: beta     = 240.       ! Constant in eq.3 (Wang et al. 2017)
! real(dp), parameter :: beta_inv = 1./beta    ! Inverse of beta constant in eq.3 (Wang et al. 2017)
! real(dp), parameter :: cstar    = 0.41       ! Unit carbon cost for maintenance of electron transport capacity (unitless) (Wang et al. 2017)
! real(dp), parameter :: A        = -3.719     ! Constant in eq.13 for water viscosity (Wang et al. 2017)
! real(dp), parameter :: B        = 580.       ! Constant in eq.13 for water viscosity (Wang et al. 2017)
! real(dp), parameter :: C        = -138.      ! Constant in eq.13 for water viscosity (Wang et al. 2017)
! real(dp), parameter :: sr2par   = 0.45       ! Conversion factor from solar irridiance (W m-2) to photosynthetically active radiation (PAR; 400-700 um) (W m-2) (Meek et al., 1984)
! real(dp), parameter :: sr2ppfd  = 2.04       ! Conversion factor from solar irradiance (W m-2) to PPFD (umol m-2 s-1) (Meek et al., 1984)
! real(dp), parameter :: coef_T   = 0.0545     ! Coefficient for temperature in eq.1 (Wang et al. 2017)
! real(dp), parameter :: coef_D   = -0.5       ! Coefficient for vapor pressure deficit in eq.1 (Wang et al. 2017)
! real(dp), parameter :: coef_z   = -0.0815    ! Coefficient for elevation in eq.1 (Wang et al. 2017)
! real(dp), parameter :: intp_C   = 1.189      ! Intercept constant in eq.1 (Wang et al. 2017)
! real(dp), parameter :: dhc      = 79430.     ! Activation energy for Kc (J mol-1) in eq.16 (Wang et al. 2017) --> from Bernacchi et al. 2001
! real(dp), parameter :: dho      = 36380.     ! Activation energy for Ko (J mol-1) in eq.16 (Wang et al. 2017) --> from Bernacchi et al. 2001
! real(dp), parameter :: kc25     = 39.97      ! k25 parameter for Kc variable (Bernacchi et al. 2001)
! real(dp), parameter :: ko25     = 27480.     ! k25 parameter for Ko variable (Bernacchi et al. 2001)
!
! ! Parameters copied rom LPJ-LMFire (gppmod and photosynthesismod)
! real(dp), parameter :: alphaa    =    0.5     ! Fraction of PAR assimilated at ecosystem level
! real(dp), parameter :: alphac3   =    0.08    ! intrinsic quantum efficiency of CO2 uptake in C3 plants
! real(dp), parameter :: alphac4   =    0.053   ! C4 intrinsic quantum efficiency
! real(dp), parameter :: bc3       =    0.015   ! leaf respiration as fraction of Vmax for C3 plants
! real(dp), parameter :: bc4       =    0.02    ! leaf respiration as fraction of vmax for C4 plants
! real(dp), parameter :: cmass     =   12.      ! atomic mass of carbon
! real(dp), parameter :: cq        =    4.6e-6  ! conversion factor for solar radiation at 550 nm from J/m2 to E/m2 (E=mol quanta)
! real(dp), parameter :: e0        =  308.56    ! parameter in Arrhenius temp response function
! real(dp), parameter :: lambdamc3 =    0.8     ! optimal (maximum?) ci:ca ratio for C3 plants
! real(dp), parameter :: lambdamc4 =    0.4     ! optimal ci:ca ratio for c4 plantsfpc_gr
! real(dp), parameter :: n0        =    7.15    ! leaf N concentration (mg/g) not involved in photosynthesis
! real(dp), parameter :: P         =    1.e5    ! atmospheric pressure (Pa)
! real(dp), parameter :: po2       =   20.9e3   ! O2 partial pressure (Pa)
! real(dp), parameter :: q10kc     =    2.1     ! q10 for temperature-sensitive parameter kc
! real(dp), parameter :: q10ko     =    1.2     ! q10 for temperature-sensitive parameter ko
! real(dp), parameter :: q10tau    =    0.57    ! q10 for temperature-sensitive parameter tau
! real(dp), parameter :: t0c3      =  250.      ! base temperature (K) in Arrhenius temperature response function for C3 plants
! real(dp), parameter :: t0c4      =  260.      ! base temperature in Arrhenius func for C4 plants
! real(dp), parameter :: tau25     = 2600.      ! value of tau at 25 deg C
! real(dp), parameter :: theta     =    0.7     ! colimitation (shape) parameter
! real(dp), parameter :: tk25      =  298.15    ! 25 deg C in Kelvin
! real(dp), parameter :: tmc3      =   45.      ! maximum temperature for C3 photosynthesis
! real(dp), parameter :: tmc4      =   55.      ! maximum temperature for C4 photosynthesis
!
! !-------------------------
! ! Pointer variables
! real(dp), pointer :: Ratm                 ! Relative atmospheric pressure to sea-level (fraction)
! real(dp), pointer :: Patm                 ! Atmospheric pressure to sea-level (Pa)
! real(dp), pointer :: vpd                  ! Average daytime saturation vapor pressure deficit (Pa)
! real(dp), pointer :: dayl                 ! Daylength (h)
! real(dp), pointer :: srad                 ! Downwelling surface shortwave radiation (kJ m-2 d-1)
! real(dp), pointer :: srad_dir             ! Direct beam downwelling shortwave raditaion (kJ m-2 d-1)
! real(dp), pointer :: srad_dif             ! Diffuse downwelling shortwave raditaion (kJ m-2 d-1)
! real(dp), pointer :: tmean                ! 24 hour mean temperature (degC)
! real(dp), pointer :: dpet                 ! Daily potential evapotranspiration (mm d-1)
! real(dp), pointer :: daet
! real(dp), pointer :: alpha
!
! real(dp), pointer :: chi                  ! Actual ratio of lead internal to external CO2 (fraction)
! real(dp), pointer :: chi0                 ! Optimal ratio of lead internal to external CO2 (fraction)
! real(dp), pointer :: gammastar            ! Photorespiratory compensation point (gammastar in Wang et al. 2017) (Pa)
! real(dp), pointer :: gpp0                 ! Gross primary productivity under non-water stressed condition (g C m-2 d-1)
! real(dp), pointer :: gpp1                 ! Gross primary productivity under actual condition (g C m-2 d-1)
! real(dp), pointer :: gpp_tot              ! Total daily gross primary productivity (g C d-1)
! real(dp), pointer :: rd                   ! Daily leaf respiration (gC m-2 d-1) >> whole day include day + night
! real(dp), pointer :: dgp                  ! Optimal daily canopy conductance (mm m-2 s-1) ??? NEED CONFIRMATION ON UNIT in LPJ-LMFire
! real(dp), pointer :: dgc                  ! Actual daily canopy conductance (mm m-2 s-1) ??? NEED CONFIRMATION ON UNIT in LPJ-LMFire
!
! real(dp), pointer :: dwscal
! real(dp), pointer :: dphen
! real(dp), pointer :: dphen_t
! real(dp), pointer :: dphen_w
!
! real(dp), pointer :: water
!
! real(dp), pointer :: fpc_grid
!
! real(dp), pointer :: elev                 ! Elevation (m)
! real(dp), pointer :: cellarea             ! Area of gridcell (m2)
! real(dp), pointer :: areafrac             ! Ground area fraction in gridcell (fraction)
!
! !-------------------------
! ! Local variables
! real(dp) :: tmean_K               ! Daily mean temperature (K)
! real(dp) :: lai                   ! Leaf area index (m2 m-2)
!
! real(dp) :: xi                    ! Constant term in Wang et al. (2017) Eq. 27
! real(dp) :: ci                    ! Optimal lead internal partial pressure of CO2 (Pa)
! real(dp) :: ca                    ! Optimal external (ambient) partial pressure of CO2 (Pa)
! real(dp) :: pi
! real(dp) :: pa
! real(dp) :: Kmm                   ! Effective Michaelis–Menten coefficient of Rubisco (Pa)
! real(dp) :: Kc                    ! Michaelis constant for CO2 in eq.16 (Wang et al. 2017)
! real(dp) :: Ko                    ! Michaelis constant for O2 in eq.16 (Wang et al. 2017)
! real(dp) :: Po                    ! O2 partial pressure in eq.17 (Pa) (Wang et al. 2017)
! real(dp) :: visco                 ! Viscosity of water relative to 25C / 298.15K (fraction)
! real(dp) :: par                   ! Photosynthetically active radiation (W m-2)
! real(dp) :: fpar                  ! Fraciton of PAR intercepted by foliage (fraction)
! real(dp) :: apar                  ! Absorbed photosynthetically active radiation (W m-2)
! real(dp) :: ppfd                  ! Photosynthetic photo flux density (mol m-2 s-1)
! real(dp) :: fppfd                 ! Fraciton of PPFD intercepted by foliage (fraction)
! real(dp) :: appfd                 ! Absorbed (by leaf) photosynthetic photo flux density (mol m-2 s-1)
! real(dp) :: chi0_term
! real(dp) :: m
!
! real(dp) :: rootprop              ! Root proportion in soil layer (fraction)
! real(dp) :: gminp                 ! Minimum canopy conductance rate (mm m-2 s-1)
!
! real(dp) :: tau
! real(dp) :: tsecs
! real(dp) :: supply                ! Supply of soil moisture (mm d-1)
! real(dp) :: demand                ! Demand for soil moisture (mm d-1)
! real(dp) :: vm0                   ! Optimal rubisco activity (gC m-2 d-1)
! real(dp) :: rd0                   ! Optimal leaf respiration rate (gC m-2 d-1)
! real(dp) :: and0                  ! Optimal daily net photosynthesis (gC m-2 d-1)
! real(dp) :: adt0                  ! Optimal total daytime net photosynthesis (gC m-2 d-1)
! real(dp) :: adtmm0                ! Optimal total daytime net photosynthesis in mm (mm m-2 d-1)
! real(dp) :: vm1                   ! Actual leaf rubisco activity (gC m-2 d-1)
! real(dp) :: adtmm1                ! Actual total daytime net photosynthesis in mm (mm m-2 d-1)
! real(dp) :: je                    ! APAR-limited photosynthesis rate (molC m-2 h-1)
! real(dp) :: jc                    ! Rubisco-activity-limited photosynthesis rate (molC m-2 h-1)
!
! !-------------------------
! ! Temperature inhibition variables
! real(dp) :: inhibx1
! real(dp) :: inhibx2
! real(dp) :: inhibx3
! real(dp) :: inhibx4
! real(dp) :: k1
! real(dp) :: k2
! real(dp) :: k3
! real(dp) :: low
! real(dp) :: high
! real(dp) :: tstress
!
! ! Rubisco capcity variables
! real(dp) :: c1
! real(dp) :: c2
! real(dp) :: b0
! real(dp) :: t0
! real(dp) :: s
! real(dp) :: sigma
!
! ! Bisection root finding variables
! real(dp) :: gpd
! real(dp) :: epsilon
! real(dp) :: x1
! real(dp) :: x2
! real(dp) :: rtbis
! real(dp) :: dx
! real(dp) :: fmid
!
! integer :: it
!
! !-------------------------
!
! Ratm     => dayvars(grid,day)%Ratm
! Patm     => dayvars(grid,day)%Patm
! vpd      => dayvars(grid,day)%vpd
! dayl     => dayvars(grid,day)%dayl
! srad     => dayvars(grid,day)%srad
! srad_dir => dayvars(grid,day)%srad_dir
! srad_dif => dayvars(grid,day)%srad_dif
! tmean    => dayvars(grid,day)%tmean
!
! elev     => topovars(grid)%elev
! cellarea => topovars(grid)%cellarea
! areafrac => topovars(grid)%areafrac
!
! chi     => vegvars(grid,day)%chi
! chi0     => vegvars(grid,day)%chi0
! gammastar      => vegvars(grid,day)%gammastar
! gpp0     => vegvars(grid,day)%gpp0
! gpp1     => vegvars(grid,day)%gpp
! gpp_tot  => vegvars(grid,day)%gpp_tot
! rd => vegvars(grid,day)%rd
!
! dgp     => vegvars(grid,day)%dgp
! dgc     => vegvars(grid,day)%dgc
! dwscal => vegvars(grid,day)%dwscal
! dphen => vegvars(grid,day)%dphen
! dphen_t => vegvars(grid,day)%dphen_t
! dphen_w => vegvars(grid,day)%dphen_w
! fpc_grid => vegvars(grid,1)%fpc_grid
!
! water => vegvars(grid,day)%water
!
! dpet => dayvars(grid,day)%dpet
! daet  => dayvars(grid,day)%daet
! alpha => dayvars(grid,day)%alpha
!
! !-------------------------
!
! lai = vegvars(grid,day)%lai
!
! if (lai == -9999.) lai = 0.
!
! tmean_K = tmean + Tfreeze
!
! if (elev < 0.) elev = 0.
! if (vpd <= 0.) vpd = 1.
! if (areafrac < 0.) areafrac = 0.
!
! dphen = 1.0
! dphen_t = 1.0
! dphen_w = 1.0
!
! !---------------------------------------------------------------------
! !---------------------------------------------------------------------
! !---------------------------------------------------------------------
! ! STEP 1: Estimation of optimal leaf CO2 partial pressure ratio from Wang et al. (2017)
! !         and daily GPP based on optimal condition
!
! !-------------------------
! ! Calculate ratio of internal to external partial pressure of CO2 in leaf from logit regression model (Eq. 1; Wang et al. 2017)
! ! T in degC; VPD in kPa; Elevation in km
!
! chi0_term = coef_T * (tmean - 25.) + coef_D * log(vpd / 1.e3) + coef_z * elev * 1.e-3 + intp_C
!
! chi0 = exp(chi0_term) / (1.0 + exp(chi0_term))
!
! ! TESTING ======= LPJ constant
! chi0 = 0.8
!
! !-------------------------
! ! Calculate photorespiratory compensation pointer (Pa)
!
! gammastar = 4.332 * Ratm * exp((37830. / 8.314) * (1./298.15  - 1./tmean_K))
!
! !-------------------------
! ! Calculate ratio of water viscosity in reference to 25degC (Eq. 13 in Wang et al. 2017)
!
! visco = exp(A + B / (C + tmean_K)) /  exp(A + B / (C + 298.15))
!
! !-------------------------
! ! Calculate photosynthetic photon flux density (PPFD) (mol m-2 s-1) for method in Wang et al. (2017)
! ! NOTE: Generalized conversion from solar irridiance used here (from Meek et al. 1984, and verified by other references; see above)
! ! Calculate PPFD absorbed by leaf based on LAI input data >> Beer's law of extinction (Eq. 1, Haxeltine & Prentice, 1996)
! ppfd = 1.e-6 * sr2ppfd * (srad_dir * 1.e3 / (dayl*3600.))      ! from micromol to mol
!
! fppfd = fpc_grid
!
! appfd = ppfd * fppfd
!
! ! appfd = ppfd - (ppfd * exp(-0.5 * lai))
!
! !-------------------------
! ! Calculate photosynthetically active radiation (PAR) (W m-2) for LPJ method in Sitch et al. (2003)
! ! Assume fPAR to be equal to fpc_grid (i.e. fraction of ground covered by foliage, thus absorbing PAR) (Leo Lai Oct 2021)
! ! alphaa = constant for fraction of PAR absorbed at ecosystem level
! par = sr2par * (srad_dir * 1.e3 / (dayl * 3600.))
!
! fpar = fpc_grid
!
! apar = par * fpar * alphaa
!
! apar = apar * daysec            ! from W m-2 to J m-2 d-1
!
!
! !-------------------------
! ! Calculate Effective Michaelis–Menten coefficient of Rubisco (Pa)
!
! Kc = kc25 * exp(dhc * (tmean_K - 298.15) / (298.15 * R * tmean_K))
! Ko = ko25 * exp(dho * (tmean_K - 298.15) / (298.15 * R * tmean_K))
!
! Po = 21000. * exp(-0.114 * elev * 1.e-3)      ! Wang et al. (2017)
! ! Po = 20947. * (Patm * 1.e-6)                ! Bernacchi et al. (2001)
!
! Kmm = Kc * (1.0 + Po / Ko)
!
! !-------------------------
! ! Calculate internal and external partial pressure of CO2 (Eq. 27 in Wang et al. 2017)
!
! xi = sqrt(beta * (Kmm + gammastar) / (1.6 * visco))
!
! ! pi = (gammastar * sqrt(vpd)) / (xi + sqrt(vpd) - xi / chi0)
! !
! ! pa = pi / chi0
!
! pa = gammastar * (1. - xi / (xi + sqrt(vpd))) / (chi0 - xi / (xi + sqrt(vpd)))
!
! pa = 28.8     ! Pa, from ghg_fix1850.nc 285 ppm
!
! pi = (xi * pa + gammastar * sqrt(vpd)) / (xi + sqrt(vpd))
!
!
! !-------------------------
! ! Calculate GPP (Eq. 2 and 3 in Wang et al. 2017)
! ! NOTE: Estimated global GPP at 123 +/- 8 Pg C per year (Beer et al., 2010)
!
! m = (pa - gammastar) / (pa + 2.0 * gammastar + 3 * gammastar * sqrt(1.6 * visco * (vpd / 1.e3) * beta_inv / (Kmm + gammastar)))
!
! gpp0 = phi0 * appfd * m * sqrt(1.0 - (cstar / m) ** 0.66)        ! g C m-2 s-1
!
! gpp0 = gpp0 * daysec                                 ! g C m-2 d-1
!
!
! !---------------------------------------------------------------------
! !---------------------------------------------------------------------
! !---------------------------------------------------------------------
! ! STEP 2: Calculate of optimal condition canopy conductance based on Haxeltine and Prentice (1996) BIOME3 paper
! !         code adpated from LPJ-LMFire gpp and photosynthesis subroutines (Leo Lai, Aug 2021)
! !
! !         Also initialize / recalculate some variables (mostly same but different parameterization) used in LPJ (Sitch et al., 2003)
!
! !-------------------------
! ! Get mole fraction of CO2 from ambient CO2 pressure calculated in previous step
! pa = pa
! pi = pi          ! In partial pressure (Pa) now after method in Wang et al. (2017)
! ca = pa / p      ! Convert to mole fraction
! ci = pi / p
!
! !-------------------------
! ! Some PFT specific parameters yet to be more precisely defined....
! inhibx1 = 2.0           ! Temperature inhibition functions
! inhibx2 = 25.0
! inhibx3 = 30.0
! inhibx4 = 55.0
!
! rootprop = 0.9          ! Root proportion
! gminp = 0.3             ! Minimum canopy conductance (mm s-1)
!
! tsecs = dayl * 3600.
!
! !-------------------------
! ! Calculate temperature inhibition function on photosynthesis
! ! High-temperature inhibition modelled conservatively as a step function prohibiting photosynthesis above 55 deg C (Table 3.7, Larcher 1983)
! ! NOTE: C3 plant assupmtion for now
!
! if (tmean < inhibx4) then
!
!    k1  = 2. * alog((1. / 0.99) - 1.) / (inhibx1 - inhibx2)
!    k2  = (inhibx1 + inhibx2) / 2.
!    low = 1./(1.+ exp(k1*(k2-tmean)))
!
!    k3   = alog(0.99 / 0.01) / (inhibx4 - inhibx3)
!    high = 1. - 0.01 * exp(k3 * (tmean - inhibx3))
!
!    tstress = low * high
!
! else    ! High temperature stress on photosynthesis, tstress = 0.
!
!    tstress = 0.
!
!    k1  = 2. * alog((1. / 0.99) - 1.) / (inhibx1 - inhibx2)
!    k2  = (inhibx1 + inhibx2) / 2.
!    low = 1./(1.+ exp(k1*(k2-tmean)))
!
!    k3   = alog(0.99 / 0.01) / (inhibx4 - inhibx3)
!    high = 1. - 0.01 * exp(k3 * (tmean - inhibx3))
!
!    tstress = low * high
!
! end if
!
! !-------------------------
! ! Recalculate parameters from LPJ paper method (Leo Lai Oct 2021)
! Ko  = ko25  *  q10ko**((tmean - 25.) / 10.)  !Michaelis constant of rubisco for O2
! Kc  = kc25  *  q10kc**((tmean - 25.) / 10.)  !Michaelis constant for CO2
! tau = tau25 * q10tau**((tmean - 25.) / 10.)  !CO2/O2 specificity ratio
!
! !-------------------------
! ! Calculation of optimal rubisco capacity, Vm (gC/m2/day) from BIOME3 model (Haxeltine & Prentice, 1996)
! c1 = tstress * alphac3 * ((pi - gammastar) / (pi + 2.0 * gammastar))    ! Eq. 4 (Haxeltine & Prentice, 1996)
!
!   ! if (tmean > tmc3) c1 = 0.
!
! c2 = (pi - gammastar) / (pi + Kc * (1.0 + Po / Ko))                     ! Eq. 6 (Haxeltine & Prentice, 1996)
!
!   b0 = bc3
!   t0 = t0c3
!
! s = (24. / dayl) * b0                                                   ! Eq. 13 (Haxeltine & Prentice, 1996)
!
! sigma = sqrt(max(0., 1. - (c2 - s)/(c2 - theta * s)))                   ! Eq. 12 (Haxeltine & Prentice, 1996)
!
! vm0 = (1. / b0) * (c1 / c2) * ((2. * theta - 1.) * s - &                ! Eq. 11 (Haxeltine & Prentice, 1996)
!      (2. * theta * s - c2) * sigma) * apar * cmass * cq
!
!
! ! Optimal daily leaf respiration, Rd, gC m-2 d-1                                ! Eq. 10 (Haxeltine & Prentice, 1996)
! rd0 = b0 * vm0
!
! ! Optimal daily net photosynthesis (at leaf level),And, gC m-2 d-1 >> using the optimal condition estimated gpp0 from (Wang et al., 2017)
! and0 = gpp0 - rd0
!
! ! Optimal total daytime net photosynthesis, Adt, gC m-2 d-1                     ! Eq. 19 (Haxeltine & Prentice, 1996)
! adt0 = and0 + (1. - dayl / 24.) * rd0
!
! ! Convert adt from gC m-2 d-1 to mm m-2 d-1 using ideal gas equation
! adtmm0 = ((adt0 / cmass) * 8.314 * (tmean + 273.3) / P) * 1000.
!
!
! !-------------------------
! ! Calculate optimal daily canopy conductance based on first estimation of photosynthesis
!
! dgp = ((1.6 * adtmm0) / (ca * (1. - chi0))) / tsecs + gminp
!
!
!
! !---------------------------------------------------------------------
! !---------------------------------------------------------------------
! !---------------------------------------------------------------------
! ! STEP 3: Calculate of temperature and water stressed GPP based on Haxeltine and Prentice (1996) BIOME3 paper
! !         code adpated from LPJ-LMFire gpp and photosynthesis subroutines (Leo Lai, Oct 2021)
!
! ! Calculate supply and demand of water on current day (mm d-1)
! supply = 5.0 * rootprop * (sum(soilvars(grid)%Tliq(1:3)) / sum(soilvars(grid)%Tsat(1:3)))
!
! ! supply = 5.0 * rootprop * sum(soilvars(grid)%Wliq(1:3)) / sum(dzmm(1:3))
!
! demand = dpet * 1.4 * (1. - exp(-dgp / 5.))
!
! ! Calculate water stress factor
! dwscal = min(1.0, supply / demand)
!
! ! Determine phenology status from water stress factor
! ! if (dwscal > 0.1) then ! .and. leafon
! !
! !   dphen_w = 1.
! !   dphen = dphen_t
! !   ! leafondays = leafondays + 1
! !
! ! else
! !
! !   dphen_w = 0.
! !   dphen = 0.
! !
! ! end if
!
! ! Second calculation of supply and demand of water based on water phenology on current day (mm d-1) (Leo Lai Oct 2021)
! supply = 5.0 * rootprop * (sum(soilvars(grid)%Tliq(1:3)) / sum(soilvars(grid)%Tsat(1:3)))
!
! ! supply = 5.0 * rootprop * sum(soilvars(grid)%Wliq(1:3)) / sum(dzmm(1:3))
!
! ! demand = dpet * 1.4 * (1. - exp(-dgp * dphen / 5.))
!
! ! Calculate water stress factor
! dwscal = min(1.0, supply / demand)
!
! daet = min(supply,demand)
!
! !------
!
! if (dpet > 0.0) then
!
!   alpha = daet / dpet
!
! else
!
!   alpha = 1.0
!
! end if
!
! !-------------------------
! ! Calculate ACTUAL canopy conductance rate from water balance conditions
! if (supply >= demand) then
!
!   dgc = dgp * dphen
!
! else if (demand > supply .and. dpet > 0.) then
!
!   dgc = max(-5. * log(1.0 - supply / (dpet * 1.4)), 0.)
!
! else
!
!   dgc = 0.0
!
! end if
!
! !-------------------------
! ! Update apar based on water based phenology (meanfpc in LPJ subroutine)
! apar = apar * dphen
!
! !---------------------------------------------------------------------
! ! Calculate the leaf CO2 partial pressure ratio in water stressed condition
! ! and then combine with temperature inhibition function for calculating actual GPP (gC m-2 d-1) and respiration (Rd, gC m-2 d-1)
!
! gpd = tsecs * (dgc - gminp)      ! Convert canopy conductance from mm m-2 s-1 to mm m-2 d-1
!
! !-------------------------
! ! Bisection root finding method to solve Eq. 2 and Eq. 18 simultaneuously (Haxeltine & Prentice, 1996)
! ! To find actual intercellular to ambient partial pressure of CO2 (chi, or lamba in LPJ)
!
! epsilon = 0.05
!
! if (gpd > 1.e-5) then     ! Significant canopy conductance
!
!   ! Initiate bracket for numerical solution
!
!   x1 = 0.02               ! Minimum bracket of the root
!   x2 = chi0               ! Maximum bracket of the root = optimal ratio
!   rtbis = x1              ! Root of the bisection
!   dx = x2 - x1
!
!   it = 0                  ! Number of tries towards solution
!
!   fmid = epsilon + 1.
!
!   do  ! Bisection root finding
!
!     it  = it + 1
!     dx  = dx * 0.5
!     chi = rtbis + dx
!
!     ! Calculate total daytime photosynthesis implied by canopy conductance from dgc calculation and
!     ! current guess for chi (xmid).
!     ! Units are mm m-2 d-1 (mm come from gpd value, mm d-1)
!     ! Reverse of Eq. 18 (Haxeltine & Prentice, 1996)
!
!     adtmm1 = (gpd * (ca * (1. - chi))) / 1.6
!
!     !Evaluate fmid at the point chi = xmid fmid will be an increasing function with xmid,
!     ! with a solution (fmid=0) between x1 and x2
!
!     fmid = adtmm0 - adtmm1
!
!     if (fmid < 0.) rtbis = chi
!
!     ! exit if solution found or > 10 iterations needed without solution
!
!     if (abs(fmid) < epsilon .or. it > 10) exit
!
!   end do ! End of bisection
!
!   !-------------------------
!   ! Calculate GPP again based on actual CO2 ratio (chi)
!   ! Find new internal leaf CO2 partial pressure assuming same external pressure from optimal calculation (Leo Lai, Aug 2021)
!
!   ! Maximum value of chi = 0.8 for C3 plant and chi = 0.4 for C4 assumed (Sitch et al., 2003; P.169)
!   chi = max(chi, 0.8)
!
!   pi = chi * pa
!
!   !-------------------------
!   ! Second calculation >>> actual rubisco capacity, Vm (gC/m2/day) from BIOME3 model (Haxeltine & Prentice, 1996)
!
!   c1 = tstress * alphac3 * ((pi - gammastar) / (pi + 2.0 * gammastar))    ! Eq. 4 (Haxeltine & Prentice, 1996)
!
!     if (tmean > tmc3) c1 = 0.
!
!   c2 = (pi - gammastar) / (pi + Kc * (1.0 + Po / Ko))                     ! Eq. 6 (Haxeltine & Prentice, 1996)
!
!     b0 = bc3
!     t0 = t0c3
!
!   s = (24. / dayl) * b0                                                   ! Eq. 13 (Haxeltine & Prentice, 1996)
!
!   sigma = sqrt(max(0., 1. - (c2 - s)/(c2 - theta * s)))                   ! Eq. 12 (Haxeltine & Prentice, 1996)
!
!   vm1 = (1. / b0) * (c1 / c2) * ((2. * theta - 1.) * s - &                 ! Eq. 11 (Haxeltine & Prentice, 1996)
!         (2. * theta * s - c2) * sigma) * apar * cmass * cq
!
!
!   ! Calculation of APAR-limited photosynthesis rate Je, molC m-2 -h;              ! Eq. 3 (Haxeltine & Prentice, 1996)
!   je = c1 * apar * cmass * cq / dayl
!
!   ! Calculation of rubisco-activity-limited photosynthesis rate Jc, molC m-2 -h   ! Eq. 5 (Haxeltine & Prentice, 1996)
!   jc = c2 * vm1 / 24.
!
!   if (je < 1.e-10 .or. jc <= 1.e-10) then
!
!     gpp1 = 0.0
!
!   else
!
!     ! Calculation of daily gross photosynthesis, Agd, gC/m2/day; Eq. 2 (Haxeltine & Prentice, 1996)
!     ! NOTE: there is an error in this equation in the above paper (missing theta in 4*theta*je*jc term) which is fixed here (LPJ-LMFire comment)
!
!     gpp1 = (je + jc - sqrt((je + jc)**2. - 4. * theta * je * jc)) / (2. * theta) * dayl
!
!   end if
!
!   ! Calculation of actual daily leaf respiration from temperature and water conditions (gC m-2 d-1)
!   rd = b0 * vm1
!
! else  !infinitesimal canopy conductance
!
!   rd   = 0.
!   gpp1 = 0.
!   chi  = 0.
!
! end if  !canopy conductance
!
! !------
!
! if (dayvars(grid,day)%dpet > 0. .and. gpp1 > 0.) then
!
!   if (year >= 1) then
!
!     gpp_tot = gpp1 * cellarea * areafrac
!
!     gpp1 = gpp1
!
!   else
!
!     gpp_tot = gpp0 * cellarea * areafrac
!
!     gpp1 = gpp0
!
!   end if
!
! end if
!
!
!
! if (areafrac == 0.) then
!   gpp0 = -9999
!   gpp_tot = 0.
! end if
!
! water = dgp !sum(soilvars(grid)%Tliq(1:6)) / sum(soilvars(grid)%Tsat(1:6))
!
! ! daet = demand
!
! ! water = soilvars(grid)%Tliq(1) / soilvars(grid)%Tsat(1)
!
! ! if (gpp1 > 0. .and. year >= 2) print *, fmid, dgp, dgc, it, chi0, gpp0, gpp_tot, rd * cellarea * areafrac, C_leaf, gtemp
!
! ! if (lprint .and. grid==gprint) write(0,*) gammastar, pi, pa, chi0, vpd, xi, gpp0, gpp1
!
!
! end subroutine gpp_new

!---------------------------------------------------------------------

end module gppmod
