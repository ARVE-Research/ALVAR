module gppmod

! Module to calculate GPP coded for ALVAR by Leo Lai (Aug 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

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
use statevarsmod,   only : ndyear,dayvars,soilvars,vegvars,topovars,gprint,lprint
use biome1mod,    only : biomevars

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

! Parameters
real(sp), parameter :: R        = 8.314      ! Universal gas constant (J mol-1 K-1)
real(sp), parameter :: phi0     = 1.02       ! Intrinsic quantum yield (g C mol-1) (Wang et al., 2017)
real(sp), parameter :: beta     = 240.       ! Constant in eq.3 (Wang et al. 2017)
real(sp), parameter :: beta_inv = 1./beta    ! Inverse of beta constant in eq.3 (Wang et al. 2017)
real(sp), parameter :: cstar    = 0.41       ! Unit carbon cost for maintenance of electron transport capacity (unitless) (Wang et al. 2017)
real(sp), parameter :: A        = -3.719     ! Constant in eq.13 for water viscosity (Wang et al. 2017)
real(sp), parameter :: B        = 580.       ! Constant in eq.13 for water viscosity (Wang et al. 2017)
real(sp), parameter :: C        = -138.      ! Constant in eq.13 for water viscosity (Wang et al. 2017)
real(sp), parameter :: sr2par   = 0.45       ! Conversion factor from solar irridiance (W m-2) to photosynthetically active radiation (PAR; 400-700 um) (W m-2) (Meek et al., 1984)
real(sp), parameter :: sr2ppfd  = 2.04       ! Conversion factor from solar irradiance (W m-2) to PPFD (umol m-2 s-1) (Meek et al., 1984)
real(sp), parameter :: coef_T   = 0.0545     ! Coefficient for temperature in eq.1 (Wang et al. 2017)
real(sp), parameter :: coef_D   = -0.5       ! Coefficient for vapor pressure deficit in eq.1 (Wang et al. 2017)
real(sp), parameter :: coef_z   = -0.0815    ! Coefficient for elevation in eq.1 (Wang et al. 2017)
real(sp), parameter :: intp_C   = 1.189      ! Intercept constant in eq.1 (Wang et al. 2017)
real(sp), parameter :: dhc      = 79430.     ! Activation energy for Kc (J mol-1) in eq.16 (Wang et al. 2017) --> from Bernacchi et al. 2001
real(sp), parameter :: dho      = 36380.     ! Activation energy for Ko (J mol-1) in eq.16 (Wang et al. 2017) --> from Bernacchi et al. 2001
real(sp), parameter :: kc25     = 39.97      ! k25 parameter for Kc variable (Bernacchi et al. 2001)
real(sp), parameter :: ko25     = 27480.     ! k25 parameter for Ko variable (Bernacchi et al. 2001)

! Parameters copied rom LPJ-LMFire (gppmod and photosynthesismod)
real(sp), parameter :: alphaa    =    0.5     ! Fraction of PAR assimilated at ecosystem level
real(sp), parameter :: alphac3   =    0.08    ! intrinsic quantum efficiency of CO2 uptake in C3 plants
real(sp), parameter :: alphac4   =    0.053   ! C4 intrinsic quantum efficiency
real(sp), parameter :: bc3       =    0.015   ! leaf respiration as fraction of Vmax for C3 plants
real(sp), parameter :: bc4       =    0.02    ! leaf respiration as fraction of vmax for C4 plants
real(sp), parameter :: cmass     =   12.      ! atomic mass of carbon
real(sp), parameter :: cq        =    4.6e-6  ! conversion factor for solar radiation at 550 nm from J/m2 to E/m2 (E=mol quanta)
real(sp), parameter :: e0        =  308.56    ! parameter in Arrhenius temp response function
real(sp), parameter :: lambdamc3 =    0.8     ! optimal (maximum?) ci:ca ratio for C3 plants
real(sp), parameter :: lambdamc4 =    0.4     ! optimal ci:ca ratio for c4 plantsfpc_gr
real(sp), parameter :: n0        =    7.15    ! leaf N concentration (mg/g) not involved in photosynthesis
real(sp), parameter :: P         =    1.e5    ! atmospheric pressure (Pa)
real(sp), parameter :: po2       =   20.9e3   ! O2 partial pressure (Pa)
real(sp), parameter :: q10kc     =    2.1     ! q10 for temperature-sensitive parameter kc
real(sp), parameter :: q10ko     =    1.2     ! q10 for temperature-sensitive parameter ko
real(sp), parameter :: q10tau    =    0.57    ! q10 for temperature-sensitive parameter tau
real(sp), parameter :: t0c3      =  250.      ! base temperature (K) in Arrhenius temperature response function for C3 plants
real(sp), parameter :: t0c4      =  260.      ! base temperature in Arrhenius func for C4 plants
real(sp), parameter :: tau25     = 2600.      ! value of tau at 25 deg C
real(sp), parameter :: theta     =    0.7     ! colimitation (shape) parameter
real(sp), parameter :: tk25      =  298.15    ! 25 deg C in Kelvin
real(sp), parameter :: tmc3      =   45.      ! maximum temperature for C3 photosynthesis
real(sp), parameter :: tmc4      =   55.      ! maximum temperature for C4 photosynthesis

!-------------------------
! Pointer variables
real(sp), pointer :: Ratm                 ! Relative atmospheric pressure to sea-level (fraction)
real(sp), pointer :: Patm                 ! Atmospheric pressure to sea-level (Pa)
real(sp), pointer :: vpd                  ! Average daytime saturation vapor pressure deficit (Pa)
real(sp), pointer :: dayl                 ! Daylength (h)
real(sp), pointer :: srad                 ! Downwelling surface shortwave radiation (kJ m-2 d-1)
real(sp), pointer :: srad_dir             ! Direct beam downwelling shortwave raditaion (kJ m-2 d-1)
real(sp), pointer :: srad_dif             ! Diffuse downwelling shortwave raditaion (kJ m-2 d-1)
real(sp), pointer :: tmean                ! 24 hour mean temperature (degC)
real(sp), pointer :: dpet                 ! Daily potential evapotranspiration (mm d-1)

real(sp), pointer :: chi                  ! Actual ratio of lead internal to external CO2 (fraction)
real(sp), pointer :: chi0                 ! Optimal ratio of lead internal to external CO2 (fraction)
real(sp), pointer :: gammastar            ! Photorespiratory compensation point (gammastar in Wang et al. 2017) (Pa)
real(sp), pointer :: gpp0                 ! Gross primary productivity under non-water stressed condition (g C m-2 d-1)
real(sp), pointer :: gpp1                 ! Gross primary productivity under actual condition (g C m-2 d-1)
real(sp), pointer :: gpp_tot              ! Total daily gross primary productivity (g C d-1)
real(sp), pointer :: rd                   ! Daily leaf respiration (gC m-2 d-1) >> whole day include day + night
real(sp), pointer :: dgp                  ! Optimal daily canopy conductance (mm m-2 s-1) ??? NEED CONFIRMATION ON UNIT in LPJ-LMFire
real(sp), pointer :: dgc                  ! Actual daily canopy conductance (mm m-2 s-1) ??? NEED CONFIRMATION ON UNIT in LPJ-LMFire
real(sp), pointer :: dwscal

real(sp), pointer :: fpc_grid

real(sp), pointer :: elev                 ! Elevation (m)
real(sp), pointer :: cellarea             ! Area of gridcell (m2)
real(sp), pointer :: areafrac             ! Ground area fraction in gridcell (fraction)

!-------------------------
! Local variables
real(sp) :: tmean_K               ! Daily mean temperature (K)
real(sp) :: lai                   ! Leaf area index (m2 m-2)

real(sp) :: ci                    ! Optimal lead internal partial pressure of CO2 (Pa)
real(sp) :: ca                    ! Optimal external (ambient) partial pressure of CO2 (Pa)
real(sp) :: Kmm                   ! Effective Michaelis–Menten coefficient of Rubisco (Pa)
real(sp) :: Kc                    ! Michaelis constant for CO2 in eq.16 (Wang et al. 2017)
real(sp) :: Ko                    ! Michaelis constant for O2 in eq.16 (Wang et al. 2017)
real(sp) :: Po                    ! O2 partial pressure in eq.17 (Pa) (Wang et al. 2017)
real(sp) :: visco                 ! Viscosity of water relative to 25C / 298.15K (fraction)
real(sp) :: par                   ! Photosynthetically active radiation (W m-2)
real(sp) :: fpar                  ! Fraciton of PAR intercepted by foliage (fraction)
real(sp) :: apar                  ! Absorbed photosynthetically active radiation (W m-2)
real(sp) :: ppfd                  ! Photosynthetic photo flux density (mol m-2 s-1)
real(sp) :: appfd                 ! Absorbed (by leaf) photosynthetic photo flux density (mol m-2 s-1)
real(sp) :: chi0_term
real(sp) :: m

real(sp) :: rootprop              ! Root proportion in soil layer (fraction)
real(sp) :: gminp                 ! Minimum canopy conductance rate (mm m-2 s-1)

real(sp) :: supply                ! Supply of soil moisture (mm d-1)
real(sp) :: demand                ! Demand for soil moisture (mm d-1)
real(sp) :: vm0                   ! Optimal rubisco activity (gC m-2 d-1)
real(sp) :: rd0                   ! Optimal leaf respiration rate (gC m-2 d-1)
real(sp) :: and0                  ! Optimal daily net photosynthesis (gC m-2 d-1)
real(sp) :: adt0                  ! Optimal total daytime net photosynthesis (gC m-2 d-1)
real(sp) :: adtmm0                ! Optimal total daytime net photosynthesis in mm (mm m-2 d-1)
real(sp) :: vm1                   ! Actual leaf rubisco activity (gC m-2 d-1)
real(sp) :: adtmm1                ! Actual total daytime net photosynthesis in mm (mm m-2 d-1)
real(sp) :: je                    ! APAR-limited photosynthesis rate (molC m-2 h-1)
real(sp) :: jc                    ! Rubisco-activity-limited photosynthesis rate (molC m-2 h-1)

!-------------------------
! Temperature inhibition variables
real(sp) :: inhibx1
real(sp) :: inhibx2
real(sp) :: inhibx3
real(sp) :: inhibx4
real(sp) :: k1
real(sp) :: k2
real(sp) :: k3
real(sp) :: low
real(sp) :: high
real(sp) :: tstress

! Rubisco capcity variables
real(sp) :: c1
real(sp) :: c2
real(sp) :: b0
real(sp) :: t0
real(sp) :: s
real(sp) :: sigma

! Bisection root finding variables
real(sp) :: gpd
real(sp) :: epsilon
real(sp) :: x1
real(sp) :: x2
real(sp) :: rtbis
real(sp) :: dx
real(sp) :: fmid

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

elev     => topovars(grid)%elev
cellarea => topovars(grid)%cellarea
areafrac => topovars(grid)%areafrac

chi     => vegvars(grid,day)%chi
chi0     => vegvars(grid,day)%chi0
gammastar      => vegvars(grid,day)%gammastar
gpp0     => vegvars(grid,day)%gpp0
gpp1     => vegvars(grid,day)%gpp
gpp_tot  => vegvars(grid,day)%gpp_tot
rd => vegvars(grid,day)%rd

dgp     => vegvars(grid,day)%dgp
dgc     => vegvars(grid,day)%dgc
dwscal => vegvars(grid,day)%dwscal
fpc_grid => vegvars(grid,1)%fpc_grid

dpet => dayvars(grid,day)%dpet

!-------------------------

lai = vegvars(grid,day)%lai

if (lai == -9999.) lai = 0.

tmean_K = tmean + Tfreeze

if (elev < 0.) elev = 0.
if (vpd <= 0.) vpd = 1.
if (areafrac < 0.) areafrac = 0.

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! STEP 1: Estimation of optimal leaf CO2 partial pressure ratio from Wang et al. (2017)
!         and daily GPP based on optimal condition

!-------------------------
! Calculate ratio of internal to external partial pressure of CO2 in leaf from logit regression model (Eq. 1; Wang et al. 2017)
! T in degC; VPD in kPa; Elevation in km

chi0_term = coef_T * (tmean - 25.) + coef_D * log(vpd / 1.e3) + coef_z * elev * 1.e-3 + intp_C

chi0 = exp(chi0_term) / (1.0 + exp(chi0_term))

!-------------------------
! Calculate photorespiratory compensation pointer (Pa)

gammastar = 4.332 * Ratm * exp((37830. / 8.314) * (1./298.15  - 1./tmean_K))

!-------------------------
! Calculate ratio of water viscosity in reference to 25degC (Eq. 13 in Wang et al. 2017)

visco = exp(A + B / (C + tmean_K)) /  exp(A + B / (C + 298.15))

!-------------------------
! Calculate photosynthetic photon flux density (PPFD) (mol m-2 s-1) for method in Wang et al. (2017)
! NOTE: Generalized conversion from solar irridiance used here (from Meek et al. 1984, and verified by other references; see above)
! Calculate PPFD absorbed by leaf based on LAI input data >> Beer's law of extinction (Eq. 1, Haxeltine & Prentice, 1996)
ppfd = 1.e-6 * sr2ppfd * (srad_dir * 1.e3 / (dayl*3600.))      ! from micromol to mol

appfd = ppfd - (ppfd * exp(-0.5 * lai))

!-------------------------
! Calculate photosynthetically active radiation (PAR) (W m-2) for LPJ method in Sitch et al. (2003)
! Assume fPAR to be equal to fpc_grid (i.e. fraction of ground covered by foliage, thus absorbing PAR) (Leo Lai Oct 2021)
! alphaa = constant for fraction of PAR absorbed at ecosystem level
par = sr2par * (srad_dir * 1.e3 / (dayl*3600.))

fpar = fpc_grid

apar = par * fpar * alphaa

apar = apar * daysec            ! from W m-2 to J m-2 d-1


!-------------------------
! Calculate Effective Michaelis–Menten coefficient of Rubisco (Pa)

Kc = kc25 * exp(dhc * (tmean_K - 298.15) / (298.15 * R * tmean_K))
Ko = ko25 * exp(dho * (tmean_K - 298.15) / (298.15 * R * tmean_K))

Po = 21000. * exp(-0.114 * elev * 1.e-3)      ! Wang et al. (2017)
! Po = 20947. * (Patm * 1.e-6)                ! Bernacchi et al. (2001)

Kmm = Kc * (1.0 + Po / Ko)

!-------------------------
! Calculate internal and external partial pressure of CO2 (Eq. 27 in Wang et al. 2017)

ci = sqrt(beta * (Kmm + gammastar) / (1.6 * visco))

ca = ci / chi0

!-------------------------
! Calculate GPP (Eq. 2 and 3 in Wang et al. 2017)
! NOTE: Estimated global GPP at 123 +/- 8 Pg C per year (Beer et al., 2010)

m = (ca - gammastar) / (ca + 2.0 * gammastar + 3 * gammastar * sqrt(1.6 * visco * (vpd / 1.e3) * beta_inv / (Kmm + gammastar)))

gpp0 = phi0 * appfd * m * sqrt(1.0 - (cstar / m) ** 0.66)        ! g C m-2 s-1

gpp0 = gpp0 * daysec                                 ! g C m-2 d-1


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! STEP 2: Calculate of temperature and water stressed GPP based on Haxeltine and Prentice (1996) BIOME3 paper
!         code adpated from LPJ-LMFire gpp and photosynthesis subroutines (Leo Lai, Aug 2021)

!-------------------------
! Some PFT specific parameters yet to be more precisely defined....
inhibx1 = 2.0           ! Temperature inhibition functions
inhibx2 = 25.0
inhibx3 = 30.0
inhibx4 = 55.0

rootprop = 0.9          ! Root proportion
gminp = 0.3             ! Minimum canopy conductance (mm s-1)

!-------------------------
! NOTE: C3 plant assupmtion for now
! Calculate temperature inhibition function on photosynthesis
! High-temperature inhibition modelled conservatively as a step function prohibiting photosynthesis above 55 deg C (Table 3.7, Larcher 1983)

if (tmean < inhibx4) then

   k1  = 2. * alog((1. / 0.99) - 1.) / (inhibx1 - inhibx2)
   k2  = (inhibx1 + inhibx2) / 2.
   low = 1./(1.+ exp(k1*(k2-tmean)))

   k3   = alog(0.99 / 0.01) / (inhibx4 - inhibx3)
   high = 1. - 0.01 * exp(k3 * (tmean - inhibx3))

   tstress = low * high

else

   tstress = 0.

end if

!-------------------------
! Calculation of optimal rubisco capacity, Vm (gC/m2/day) from BIOME3 model (Haxeltine & Prentice, 1996)
c1 = tstress * alphac3 * ((ci - gammastar) / (ci + 2.0 * gammastar))    ! Eq. 4 (Haxeltine & Prentice, 1996)

  if (tmean > tmc3) c1 = 0.

c2 = (ci - gammastar) / (ci + Kc * (1.0 + Po / Ko))                     ! Eq. 6 (Haxeltine & Prentice, 1996)

  b0 = bc3
  t0 = t0c3

s = (24. / dayl) * b0                                                   ! Eq. 13 (Haxeltine & Prentice, 1996)

sigma = sqrt(max(0., 1. - (c2 - s)/(c2 - theta * s)))                   ! Eq. 12 (Haxeltine & Prentice, 1996)

vm0 = (1. / b0) * (c1 / c2) * ((2. * theta - 1.) * s - &                ! Eq. 11 (Haxeltine & Prentice, 1996)
     (2. * theta * s - c2) * sigma) * apar * cmass * cq


! Optimal daily leaf respiration, Rd, gC m-2 d-1                                ! Eq. 10 (Haxeltine & Prentice, 1996)
rd0 = b0 * vm0

! Optimal daily net photosynthesis (at leaf level),And, gC m-2 d-1 >> using the optimal condition estimated gpp0 from (Wang et al., 2017)
and0 = gpp0 - rd0

! Optimal total daytime net photosynthesis, Adt, gC m-2 d-1                     ! Eq. 19 (Haxeltine & Prentice, 1996)
adt0 = and0 + (1. - dayl / 24.) * rd0

! Convert adt from gC m-2 d-1 to mm m-2 d-1 using ideal gas equation
adtmm0 = ((adt0 / cmass) * 8.314 * (tmean + 273.3) / P) * 1000.


!-------------------------
! Calculate optimal daily canopy conductance based on first estimation of photosynthesis

dgp = (((1.6 * adtmm0) / ((ca / P) * (1. - chi0)))) / (dayl * 3600.) + gminp

! Calculate supply and demand of water on current day (mm d-1)
supply = 5.0 * rootprop * sum(soilvars(grid)%Tliq(1:6)) / 6.

demand = dpet * 1.4 * (1 - exp(-dgp / 5.))

! Calculate water stress factor
dwscal = min(1.0, supply / demand)

!-------------------------
! Determine the actual canopy conductance rate
if (supply >= demand) then

  dgc = dgp

else if (demand > supply .and. dpet > 0.) then

  dgc = max(-5. * log(1.0 - supply / (dpet * 1.4)), 0.)

else

  dgc = 0.0

end if


!---------------------------------------------------------------------
! Calculate the leaf CO2 partial pressure ratio in water stressed condition
! and then combine with temperature inhibition function for calculating actual GPP (gC m-2 d-1) and respiration (Rd, gC m-2 d-1)

gpd = (dayl * 3600.) * (dgc - gminp)      ! Convert from mm m-2 s-1 to mm m-2 d-1

epsilon = 0.05

!-------------------------
! Bisection root finding method to solve Eq. 2 and Eq. 18 simultaneuously (Haxeltine & Prentice, 1996)

if (gpd > 1.e-5) then     ! Significant canopy conductance

  ! Initiate bracket for numerical solution

  x1 = 0.02               ! Minimum bracket of the root
  x2 = chi0               ! Maximum bracket of the root = optimal ratio
  rtbis = x1              ! Root of the bisection
  dx = x2 - x1

  it = 0                  ! Number of tries towards solution

  fmid = epsilon + 1.

  do  ! Bisection root finding

    it  = it + 1
    dx  = dx * 0.5
    chi = rtbis + dx

    ! Calculate total daytime photosynthesis implied by canopy conductance from dgc calculation and
    ! current guess for chi (xmid).
    ! Units are mm m-2 d-1 (mm come from gpd value, mm d-1)
    ! Reverse of Eq. 18 (Haxeltine & Prentice, 1996)

    adtmm1 = (gpd * ((ca / 101325.) * (1. - chi))) / 1.6

    !Evaluate fmid at the point chi = xmid fmid will be an increasing function with xmid,
    ! with a solution (fmid=0) between x1 and x2

    fmid = adtmm0 - adtmm1

    if (fmid < 0.) rtbis = chi

    ! exit if solution found or > 10 iterations needed without solution

    if (abs(fmid) < epsilon .or. it > 10) exit

  end do ! End of bisection

  !-------------------------
  ! Calculate GPP again based on actual CO2 ratio (chi)
  ! Find new internal leaf CO2 partial pressure assuming same external pressure from optimal calculation (Leo Lai, Aug 2021)

  ci = chi * ca

  !-------------------------
  ! Second calculation >>> actual rubisco capacity, Vm (gC/m2/day) from BIOME3 model (Haxeltine & Prentice, 1996)

  c1 = tstress * alphac3 * ((ci - gammastar) / (ci + 2.0 * gammastar))    ! Eq. 4 (Haxeltine & Prentice, 1996)

    if (tmean > tmc3) c1 = 0.

  c2 = (ci - gammastar) / (ci + Kc * (1.0 + Po / Ko))                     ! Eq. 6 (Haxeltine & Prentice, 1996)

    b0 = bc3
    t0 = t0c3

  s = (24. / dayl) * b0                                                   ! Eq. 13 (Haxeltine & Prentice, 1996)

  sigma = sqrt(max(0., 1. - (c2 - s)/(c2 - theta * s)))                   ! Eq. 12 (Haxeltine & Prentice, 1996)

  vm1 = (1. / b0) * (c1 / c2) * ((2. * theta - 1.) * s - &                 ! Eq. 11 (Haxeltine & Prentice, 1996)
        (2. * theta * s - c2) * sigma) * apar * cmass * cq


  ! Calculation of APAR-limited photosynthesis rate Je, molC m-2 -h;              ! Eq. 3 (Haxeltine & Prentice, 1996)
  je = c1 * apar * cmass * cq / dayl

  ! Calculation of rubisco-activity-limited photosynthesis rate Jc, molC m-2 -h   ! Eq. 5 (Haxeltine & Prentice, 1996)
  jc = c2 * vm1 / 24.

  if (je < 1.e-10 .or. jc <= 1.e-10) then

    gpp1 = 0.0

  else

    ! Calculation of daily gross photosynthesis, Agd, gC/m2/day; Eq. 2 (Haxeltine & Prentice, 1996)
    ! NOTE: there is an error in this equation in the above paper (missing theta in 4*theta*je*jc term) which is fixed here (LPJ-LMFire comment)

    gpp1 = (je + jc - sqrt((je + jc)**2. - 4. * theta * je * jc)) / (2. * theta) * dayl

  end if

  ! Calculation of actual daily leaf respiration from temperature and water conditions (gC m-2 d-1)
  rd = b0 * vm1

else  !infinitesimal canopy conductance

  rd   = 0.
  gpp1 = 0.
  chi  = 0.

end if  !canopy conductance

!------

if (dayvars(grid,day)%dpet > 0. .and. gpp1 > 0.) then

  if (year >= 1) then

    gpp1 = gpp1
    gpp_tot = gpp1 * cellarea * areafrac

  else

    gpp1 = gpp0
    gpp_tot = gpp0 * cellarea * areafrac

  end if

end if



if (areafrac == 0.) then
  gpp0 = -9999
  gpp_tot = 0.
end if


! if (gpp1 > 0. .and. year >= 2) print *, fmid, dgp, dgc, it, chi0, gpp0, gpp_tot, rd * cellarea * areafrac, C_leaf, gtemp


end subroutine gpp

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
! real(sp), parameter, dimension(17) :: sla0 = [0.012    &       ! 1 = Tropical evergreen / Broadlead evergreen tree (tropical)
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
! real(sp), parameter, dimension(17) :: m    = [0.0015    &       ! 1 = Tropical evergreen / Broadlead evergreen tree (tropical)
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
! real(sp), pointer :: C_leaf
! real(sp), pointer :: sla
! real(sp), pointer :: lai
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
! real(sp), parameter, dimension(17) :: sla0 = [0.012    &       ! 1 = Tropical evergreen / Broadlead evergreen tree (tropical)
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
! real(sp), parameter, dimension(17) :: m    = [0.0015    &       ! 1 = Tropical evergreen / Broadlead evergreen tree (tropical)
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
! real(sp), pointer :: C_leaf
! real(sp), pointer :: sla
! real(sp), pointer :: lai0
!
! real(sp) :: lai_sun
! real(sp) :: lai_sha
!
! real(sp), parameter :: kb = 0.5       ! Extinction factors, should be function of zenith angle and leaf angle, but constant FOR NOW (Dai et al., 2004)
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

!---------------------------------------------------------------------

end module gppmod
