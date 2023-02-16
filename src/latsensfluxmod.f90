module latsensfluxmod

implicit none

contains

!-----------------------------------------------------------------------
subroutine latsensflux(i,snl,dpet,tday,tnight,tdew,cldf,wind,rhum,Patm,Patm30,srad,fsnow,zsno,raw_b,rah_b,theta_fc,&
                       Tsoil,Tsoiln,Psat,Tsat,dz,Bexp,Wliq,Wice,&
                       Lwg_b,dLgdt_b,Hg,dHgdT,Hg_b,dHgdt_b,Eg,dEgdT,Eg_b,dEgdt_b,Beta_soil,dqsath_dT,dqgdTg,qg,qs,lam)

! Calculate the 'equilibrium' potential evapotranspiration and sensible
! heat flux for bare soil, ground under canopy, or vegetated surfaces

!NOTE: At present, ARVE uses the CLM method for calculation of sensible
!heat fluxes but uses a variant of the Prescott equation (Prentice et al.
!2003) for latent heat flux. This approach for the latent heat fluxes was chosen
!as we do not presently take in the relative humidity of the location. One
!weakness of this approach is that dew drop is not possible. The code is set up
!to allow an easy change for inclusion of a more robust approach to latent heat fluxes.
!Further improvements to ARVE should include RH in the weather generator . JM 24.05.2011

use parametersmod, only : i4,sp,dp,pi,Tfreeze
use pftparmod,     only : npft
use statevarsmod,  only : ns,nl

implicit none

! real(dp), pointer, dimension(:,:) :: temp_veg           !timestep average vegetation temperature (K)
! real(dp), pointer, dimension(:) :: pet                  !potential evapotranspiration for canopy (mm s-1)
! real(dp), pointer, dimension(:) :: red_pet              !potential evapotranspiration for canopy reduced by canopy evaporation (mm s-1)
integer, intent(in) :: i
integer, intent(in) :: snl                                 !soil layers
real(sp), intent(in) :: dpet                               !potential evapotranspiration for bare ground with no veg (mm s-1)
real(sp), intent(in) :: tday
real(sp), intent(in) :: tnight
real(sp), intent(in) :: tdew
real(sp), intent(in) :: cldf
real(sp), intent(in) :: wind                               ! Windspeed (m s-1)
real(sp), intent(in) :: rhum                               ! Relative humidity (%)
! real(dp), pointer, dimension(:) :: underpet             !potential evapotranspiration for ground beneath canopy (mm s-1)
! real(dp), pointer, dimension(:) :: betaT                !soil moisture function limiting transpiration (fraction)
! real(dp), pointer, dimension(:) :: lai_sno              !leaf area index (with snow burial)
! real(dp), pointer, dimension(:) :: stem_sno             !stem "                "        "
! real(dp), pointer, dimension(:) :: lai_sun              !leaf area index (sunlight)
! real(dp), pointer, dimension(:) :: lai_sha              !leaf are index (shaded)
real(sp), intent(in) :: Patm                               ! atmospheric pressure (Pa)
real(sp), intent(in) :: Patm30                               ! atmospheric pressure (Pa)
real(sp), intent(in) :: srad
real(sp), intent(in) :: fsnow        !fraction of the gridcell covered by snow (fraction)
real(sp), intent(in) :: zsno                               !snow depth (m)
real(sp), intent(in) :: raw_b                              !resistance to latent heat transfer (bare ground and atm)(s m-1)
! real(dp), pointer, dimension(:) :: raw                  !        "        "        "        (canopy and atm)(s m-1)
! real(dp), pointer, dimension(:) :: rawp                 !        "        "        "        (canopy and ground)(s m-1)
real(sp), intent(in) :: rah_b                              !resistance to sensible heat transfer (bare ground and atm)(s m-1)
! real(sp), dimension(:) :: rah                  !        "        "        "        (canopy and atm)(s m-1)
! real(sp), dimension(:) :: rahp                 !        "        "        "        (canopy and ground)(s m-1)
! real(sp), dimension(:) :: rdbp                 !fraction of potential evapotranspiration
! real(sp), dimension(:) :: rb                   !leaf boundary layer resistance (s m-1)
! real(sp), dimension(:) :: Wcan                 !water in the canopy (mm)
! real(sp), dimension(:) :: f_wet                !fraction of the canopy that is wet
! real(sp), dimension(:) :: f_dry                !fraction of the canopy that is dry
! real(sp), dimension(:) :: gt_sun               !canopy conductance (from last time step) (mm/s) sunlight
! real(sp), dimension(:) :: gt_sha               !canopy conductance (from last time step) (mm/s) shaded
! real(dp), pointer :: theta_fc                           !the volumetric water content at field capacity (fraction)
! real(dp), pointer :: ustar_b                            !friction velocity estimated from scalar wind speed (m s-1)

! real(dp), pointer, dimension(:) :: rdp_dry              !term for contributions by wet leaves and transpiration, limited by avail H2O and pot. evapo.
real(sp), dimension(nl) , intent(in)   :: theta_fc     ! Soil water volumetric content at field capacity (Psi = -33 kPa)   (fraction / m3 m-3)
real(sp), dimension(ns:nl), intent(in) :: Tsoil                !soil temperature (K)
real(sp), dimension(ns:nl), intent(in) :: Tsoiln
real(sp), dimension(ns:nl), intent(in) :: Psat                 !
real(sp), dimension(ns:nl), intent(in) :: Tsat                 !
real(sp), dimension(ns:nl), intent(in) :: dz                   !
real(sp), dimension(ns:nl), intent(in) :: Bexp                 !
real(sp), dimension(ns:nl), intent(in) :: Wliq                 !
real(sp), dimension(ns:nl), intent(in) :: Wice                 !
! real(dp), pointer, dimension(:) :: can_evap             !canopy evaporation (mm s-1)
real(sp), intent(inout) :: Lwg_b
real(sp), intent(inout) :: dLgdT_b
real(sp), intent(inout) :: Hg
real(sp), intent(inout) :: dHgdT
real(sp), intent(inout) :: Hg_b
real(sp), intent(inout) :: dHgdt_b
real(sp), intent(inout) :: Eg
real(sp), intent(inout) :: dEgdT
real(sp), intent(inout) :: Eg_b
real(sp), intent(inout) :: dEgdt_b
real(sp), intent(inout) :: Beta_soil                         ! CLM 4.0 Eqn 5.68 represent molecular diffusion processes in evaporation from soil
real(sp), intent(inout) :: dqsath_dT                         ! Saturated water vapour specific humidity (CLM 4 Eqn 5.143) derivative w.r.t grnd temp.
real(sp), intent(inout) :: dqgdTg                            ! Derivative of the specific humidity of the soil surface (CLM 4.0 Eqn 5.74)
real(sp), intent(inout) :: qg                                ! Ground specific humidity
real(sp), intent(inout) :: qs                                ! Canopy specific humidity
real(sp), intent(inout) :: lam

!parameter
real(sp), parameter :: grav  = 9.80616       !Gravitational acceleration  (m s-2)
real(sp), parameter :: Pstd  = 101325.       !Standard pressure           (Pa)
real(sp), parameter :: sigsb = 5.67e-8          !Stefan-Boltzmann constant   (W m-2 K-4)
real(sp), parameter :: kbz   = 1.38065e-23      !Boltzmann constant          (J K-1 molecule-1)
real(sp), parameter :: avo   = 6.0221415e26     !Avogadro's number           (molecule kmol-1)
real(sp), parameter :: Rgas  = avo * kbz        !Universal gas constant      (J K-1 kmol-1)
real(sp), parameter :: MWda  = 28.966       !Molecular weight of dry air (kg kmol-1)
real(sp), parameter :: Rda   = Rgas / MWda      !Dry air gas constant        (J K-1 kg-1)
real(sp), parameter :: MWwv  = 18.016        !Mol. weight of water vapor  (kg kmol-1)
real(sp), parameter :: Rwv   = Rgas / MWwv      !Water vapor gas constant    (J K-1 kg-1)
real(sp), parameter :: Cair  = 1.25e3        ! Heat capacity of dry air (J kg-1 K-1) (CLM parameterization)
real(sp), parameter :: pliq  = 1000.         !density of water (kg m-3)
real(sp), parameter :: pice  =  917.         !density of ice (kg m-3)
real(sp), parameter :: lvap  = 2.501e6          !Latent heat of vaporization (J kg-1)
real(sp), parameter :: Lf    = 3.337e5          !Water latent heat of fusion (J kg-1)
real(sp), parameter :: lsub  = lvap + Lf        !Latent heat of sublimation  (J kg-1)
real(sp), parameter :: asnow = 1.0                   !alpha factor for soil humidity calculation (snow)
real(sp), parameter :: s_psi_max = -1.e8                !limit for surface soil matric potential (mm)
real(sp), parameter :: dz_litter = 0.05              !assumed typical litter depth (m)
real(sp), parameter :: leaf_lit = 0.5                !assumed effective litter area index (m2 m-2)
real(sp), parameter :: minplant  =  0.05  !threshold plant and stem area for energy calculations
real(sp), parameter :: Esoil  = 0.96         !thermal emissivity of bare soil (fraction)
real(sp), parameter :: Ewater = 0.96         !thermal emissivity of water (fraction)
real(sp), parameter :: Esnow  = 0.97         !thermal emissivity of snow (fraction)
real(sp), parameter :: la = 10.77             !coefficients for the correction of downwelling longwave dependent on cld cover
real(sp), parameter :: lb =  2.34             !
real(sp), parameter :: lc =-18.44             !

!local variables
integer :: pft
real(sp) :: dt
real(sp) :: gpet
real(sp) :: TairK
real(sp) :: TdewK
real(sp) :: eatm
real(sp) :: ratm
real(sp) :: theta
real(sp) :: ustar_b
real(sp) :: gamma                                       !psychrometer constant (Pa K-1)
real(sp) :: ss                                          !rate of increase of saturated vapor pressure with temperature (Pa K-1)
real(sp) :: gamma_soil                                  !psychrometer constant (Pa K-1) for soil temperature
real(sp) :: ss_soil                                     !rate of increase of saturated vapor pressure with (soil) temperature (Pa K-1)
real(sp) :: qatm                                        !atmos. specific humidity (kg kg-1)
real(sp) :: qsath                                       !saturation specific humidity (kg kg-1)
real(sp) :: alpha                                       !alpha factor for soil humidity calculation (combined)
real(sp) :: asoil                                       !alpha factor for soil humidity calculation (soil)
real(sp) :: r_stomata_sun                               !sunlight stomatal resistance (s/mm)
real(sp) :: r_stomata_sha                               !shaded stomatal resistance (s/mm)
real(sp) :: surf_wet                                    !surface soil wetness (fraction)
real(sp) :: sPsi                                        !surface soil water matric potential (mm)
real(sp) :: theta_1                                     !soil surface volumetric water content
real(sp) :: prop                                        !theta_1 / theta_fc
real(sp) :: esat_val                                    !value of esat from called function
real(sp) :: Evpot                                       ! potential evaporation from wet foliage per unit wetted area
real(sp) :: cah                                         !sensible heat conductances (canopy air to atm)
real(sp) :: cgh                                         !        "                (ground to canopy)
real(sp) :: cvh                                         !        "                (leaf to canopy)
real(sp) :: caw                                         !latent heat conductances (canopy air to atm)
real(sp) :: cgw                                         !        "                (ground to canopy)
real(sp) :: cvw                                         !        "   real(sp), intent(in) :: Patm30             (leaf to canopy)
real(sp) :: Edem_canopy                                 !demand for evaporation from the wetted canopy
real(sp) :: canopy_avail                                !water in the canopy that is available for evaporation
real(sp) :: fsno_filt                                   !snow cover of litter
real(sp) :: r_litter                                    !litter resistance (s m-1)
real(sp) :: litt_eff                                    !find effective litter area index (m2 m-2) (area not covered by snow)

real(sp) :: D
real(sp) :: Latm
real(sp) :: Egrnd

real(sp) :: Swg_b

!-----------------------
! From calcshortstep.f90 in ARVE

if (Wliq(snl+1) == 0. .and. Wice(snl+1) > 0.) then
  lam = lsub
else
  lam = lvap
end if

if (i == 1) then
  TairK = tday + Tfreeze
  TdewK = tdew + Tfreeze
else
  TairK = tnight + Tfreeze
  TdewK = tdew + Tfreeze
end if

! Atmospheric potential temperature (K)
theta = TairK * ((Patm/Patm30)**(Rda / Cair))        !eq. 12.1 p. 161

! Atmospheric vapor pressure (Pa)
eatm = rhum / 100. * esat(TairK)  !CLM p.162

! Density of moist air (kg m-3)
ratm = (Patm - 0.378 * eatm) / (Rda * TairK)  !eqn 5.8

!-----------------------
!-----------------------
!-----------------------
!-----------------------
!-----------------------
! Routine for radflux copied from radflux for now

D = TdewK - TairK
Latm = sigsb * (TairK + la * cldf**2 + lb * cldf + lc + 0.84 * (D + 4.01))**4

! calc ground emissivity
Egrnd = (Esnow * fsnow) + (Esoil * (1. - fsnow))

!Longwave fluxes on bare ground
Lwg_b = Egrnd * sigsb * Tsoiln(snl+1)**4 -    &                                    !upwelling longwave
       Egrnd * Latm + &                                                                 !downwelling
       4. * Egrnd * sigsb * Tsoiln(snl+1)**3 * (Tsoil(snl+1) - Tsoiln(snl+1))         !rate of change

dLgdT_b = -4. * Egrnd * sigsb * Tsoiln(snl+1)**3


!-----------------------
!-----------------------
!-----------------------
!-----------------------
!-----------------------

!-----------------------
! Calculate the latent heat flux for bare ground first
! Initial pre-calcs.
esat_val = esat(Tsoil(snl+1))

qsath = 0.622 * esat_val/ (Patm - 0.378 * esat_val)  !CLM 4.0 eqn 5.142

dqsath_dT = 0.622 * Patm * desdT(Tsoil(snl+1)) / ((Patm - 0.378 * esat_val)*(Patm - 0.378 * esat_val))  !CLM 4.0 eqn 5.143

! r_litter is the litter resistance, introduced in Sakaguchi and Zeng 2009.
! find snow cover of litter
fsno_filt = zsno / dz_litter

! find effective litter area index (m2 m-2) (area not covered by snow)
litt_eff = leaf_lit * (1.0 - min(fsno_filt, 1.0))

! Litter resistance (s m-1) CLM 4.0 Eqn 5.106
ustar_b  = 0.14 * wind
r_litter = 1.0 / (0.004 * ustar_b) * (1.0 - exp(-litt_eff))


!-----------------------
! Find the surface wetness and surface soil water matric potential (CLM 5.65-5.66)
surf_wet = 1.0 / (dz(1) * Tsat(1)) * ((Wliq(1) / pliq) + (Wice(1) / pice))
surf_wet = min(surf_wet, 1.0)
surf_wet = max(surf_wet, 0.01)
sPsi = max(Psat(1) * surf_wet**(-Bexp(1)), s_psi_max)

asoil = exp((sPsi * grav) / (1.e3 * Rwv * Tsoil(1)))  !CLM 5.64
alpha = asoil * (1.0 - fsnow) + asnow * fsnow   !CLM 5.63

qg = qsath * alpha

! eatm (atmospheric vapour pressure) is precalculated in calcshortstep
qatm = (0.622 * eatm) / (Patm - 0.378 * eatm)  !CLM

! add limiting condition to prevent large increase (decreases) in qg for
! small increases (decreases) in soil moisture in very dry soils (CLM 3 tech note pg .58)
if (qsath > qatm .and. qatm > qg) then
  qg = qatm
end if


!-----------------------
! CLM 4.0 added in a new Beta function to represent the molecular diffusion process
! from the soil pore to the surface within the very dry part of the soil (ref: Sakaguchi
! and Zeng 2009). This will act to limit evaporation.

! soil surface volumetric water content (CLM 4 eq. 5.69)
theta_1 = min(1.0, 1.0 / dz(1) * ((Wliq(1) / pliq) + (Wice(1) / pice)))
theta_1 = max(theta_1, 0.01)

! the volumetric water content at field capacity (theta_fc) is derived by assuming a hydraulic
! conductivity of 0.1mm/day and inverting the hydraulic conductivity function (CLM 4 eq. 5.70)
! this is precalculated in soilinit

if (theta_1 >= theta_fc(1) .or. (qatm - qg) > 0.) then
  Beta_soil = 1.0
else if (theta_1 < theta_fc(1)) then
  prop = min(theta_1 / theta_fc(1), 1.0)
  prop = max(prop, 0.01)
  Beta_soil = 0.25 * (1.0 - fsnow) * (1.0 - cos(pi * prop))**2 + fsnow
end if


!-----------------------
! Find the latent heat flux from bare soil (gpet)
!----------Prescott
!psychrometer constant, weakly temperature dependent (uses degrees C, but units of gamma are Pa K-1)
gamma_soil = 65.05 + Tsoil(snl+1) * 0.064

! Find the rate of increase of saturated vapor pressure with temperature (Pa K-1)
ss_soil = desdT(Tsoil(snl+1))

gpet = max((ss_soil / (ss_soil + gamma_soil)) * (Swg_b - Lwg_b) / lvap, 0.)  !(mm s-1)

!----------
! Not used at present.
! gpet = -ratm * Beta_soil * (qatm - qg) / raw_b  !CLM 4.0 Eqn 5.62

!-----------------------

Eg_b = gpet

dqgdTg = alpha * dqsath_dT  !CLM 4.0 Eqn 5.74

dEgdt_b = Beta_soil * ratm * dqgdTg / raw_b   !CLM 4.0 Eqn 5.73


!-----------------------
! Calculate sensible heat flux between the air and soil surface (CLM eqn 5.60, CLM 4.0 Eqn 5.61)
Hg_b = -ratm * Cair * (theta - Tsoil(snl+1)) / rah_b

dHgdt_b = ratm * Cair / rah_b

! print *, i,srad-Lwg_b-Hg_b-lam*Eg_b,srad, -Lwg_b,-Hg_b,-lam*Eg_b,-ratm,Cair,theta,Tsoil(snl+1)-273.15,&
!           tday-Tsoil(snl+1)+273.15,tnight-Tsoil(snl+1)+273.15,rah_b

! ! Now find the latent and sensible heat fluxes for vegetation and the ground below vegetation
! do pft = 1,npft
!     if (lai_sno(pft) + stem_sno(pft) == 0._dp) then !no vegetation
!
!             !no canopy hence 0
!             pet(pft) = 0._dp
!             underpet(pft) = 0._dp
!             red_pet(pft) = 0._dp
!             can_evap(pft) = 0._dp
!             Hg_v(pft) = 0._dp
!             dEgdt_v(pft) = 0._dp
!             dHgdt_v(pft) = 0._dp
!             dLgdt_v(pft) = 0._dp
!
!             cycle
!
!     else if (lai_sno(pft) + stem_sno(pft) > minplant) then !vegetated surfaces
!
!         ! Calculate vegetation potential evapotranspiration using the LPJ relation
!
!             !----------Prescott
!             !psychrometer constant, weakly temperature dependent (uses degrees C, but units of gamma are Pa K-1)
!             gamma = 65.05_dp + Tair * 0.064_dp
!
!             !find the rate of increase of saturated vapor pressure with temperature (Pa K-1)
!             ss = desdT(TairK)
!
!             ! variant of the Prescott equation (Prentice et al. 1993)
!             pet(pft) = max((ss / (ss + gamma)) * (Swv(pft) - Lwv(pft)) / lvap, 0._dp)  !(mm s-1)
!             !----------
!
!             ! Reduce the PET by the amount of water that can be evaporated from the wetted portion of the canopy.
!                      if (pet(pft) > 0._dp .and. Wcan(pft) > 0._dp) then
!
!                             ! The water flux from wetted leaves and stems can be found from the amount of water in the canopy
!                              canopy_avail = Wcan(pft) / dt !(mm s-1)
!
!                              ! The PET is decreased by the fraction of leaf and stem area that is wet and non-transpiring
!                              red_pet(pft) = pet(pft) * f_dry(pft)
!
!                              ! The demand for evaporation from the wetted canopy is then
!                              Edem_canopy = pet(pft) * f_wet(pft)
!
!                              if (canopy_avail > Edem_canopy) then
!
!                                  can_evap(pft) = Edem_canopy !(mm s-1)
!
!                              else ! the evaporative demand is more than what is available in the canopy
!
!                                  can_evap(pft) = canopy_avail   !(mm s-1)
!
!                              ! the excess demand goes back to red_pet
!                              red_pet(pft) = red_pet(pft) + Edem_canopy - canopy_avail
!
!                              end if
!
!                        else if (pet(pft) > 0._dp .and. Wcan(pft) == 0._dp) then !evaporation is not reduced and canopy is dry
!
!                              can_evap(pft) = 0._dp
!                              red_pet(pft) = pet(pft)
!
!                    else if (pet(pft) <= 0._dp) then
!
!                   ! Condensing conditions can not exist as we have no way of knowing the amount of water to condense
!
!                         can_evap(pft) = 0._dp
!                         red_pet(pft) = 0._dp
!
!                        end if !condensing/evap
!
!         ! Calculate the latent and sensible heat flux below vegetation
!
!              ! use eatm, qatm,qsath from the soil temp as we do not uniquely calc the veg temp
!
!             ! potential evaporation from wet foliage per unit wetted area
!             ! this uses the last timestep's qs value.qs is canopy specific humidity. Somewhat
!             ! meaningless at present as we do not track canopy humidity...
!             Evpot = - ratm * (qs - qsath) / rb(pft)                !CLM 4.0 Eqn 5.95
!
!             !calc fraction of potential evapotranspiration !CLM 4.0 Eqn 5.94
!             if (Evpot > 0._dp .and. betaT(pft) > 1.e-10) then
!
!                  !NOTE: this is using the stomatal resistance from the time step before.
!                  !FLAG CHECK if this is ok! JM 21.12.2010
!                 if (gt_sun(pft) > 0. .and. gt_sha(pft) > 0. .and. lai_sno(pft) > 0._dp) then
!                  r_stomata_sun = 1._dp / gt_sun(pft)
!                  r_stomata_sha = 1._dp / gt_sha(pft)
!
!                 rdp_dry(pft) = (f_dry(pft) * raw_b / lai_sno(pft)) * (lai_sun(pft) / (raw_b + r_stomata_sun) &
!                                                 + lai_sha(pft) / (raw_b + r_stomata_sha))
!
!                 rdbp(pft) = min(f_wet(pft) + rdp_dry(pft), (Evpot * rdp_dry(pft) + Wcan(pft) / dt) / Evpot)
!
!                 else
!                    rdbp(pft) = f_wet(pft)
!                 end if
!
!             else if (Evpot > 0._dp .and. betaT(pft) <= 1.e-10) then
!                rdbp(pft) = min(f_wet(pft),(Evpot * rdp_dry(pft) + Wcan(pft) / dt) / Evpot)
!             else
!                rdbp(pft) = 1._dp
!             end if
!
!             !calculate the potential evapotranspiration for that skin temperature of the ground below the vegetation
!
!             !calculate the conductances to latent heat transfer
!             caw = 1._dp / raw(pft)
!             cgw = Beta_soil / (rawp(pft) + r_litter)        !CLM 4.0 changes this to add in rlitter (see eqn. 5.93)
!             cvw = (lai_sno(pft) + stem_sno(pft)) / raw_b * rdbp(pft)
!
!             !----------Prescott
!             underpet(pft) = max((ss_soil / (ss_soil + gamma_soil)) * (Swg_v(pft) - Lwg_v(pft)) / lvap, 0._dp)  !(mm s-1)
!             !----------
!
!             !NOT presently used---
!             !canopy specific humidity (CLM 4.0 Eqn 5.90)
!             !qs = caw * qatm + cgw * qg + cvw * qsath / (caw + cvw + cgw)
!
!             !PET under the canopy
!             !underpet(pft) = -ratm * (caw * qatm + cvw * qsath - (caw + cvw) * qg) * cgw / (caw + cvw + cgw)  !CLM 5.97
!             !---
!
!             Eg_v(pft) = underpet(pft)
!
!            !find the sensible heat conductances
!                        !from the canopy air to the atmosphere
!                        cah = 1. / rah(pft)    !eq. 5.76
!
!                        !from the ground to the canopy
!                        cgh = 1. / rahp(pft)
!
!                        !from the leaf to the canopy
!                        cvh = (lai_sno(pft) + stem_sno(pft)) / rb(pft)
!
!               !calculate the sensible heat flux from the ground
!               Hg_v(pft) = -ratm * Cair * (cah * theta + cvh * temp_veg(2,pft) - &
!                     (cah + cvh) * Tsoil(snl+1)) * (cgh / (cah + cvh + cgh)) !eq. 5.80
!
!               dEgdt_v(pft) = Beta_soil * ratm / (rawp(pft) + r_litter) * (caw + cvw) / (caw + cvw + cgw) * dqgdTg  !CLM 4.0 Eqn 5.111
!
!               dHgdt_v(pft) = ratm * Cair / rahp(pft) * (cah + cvh) / (cah + cvh + cgh)  ! CLM 4.0 Eqn 5.110
!
!               ! this presently uses the dqsath_dT value for the soil, not veg temp.
!               dLgdt_v(pft) = lam * ratm * (caw + cgw) * cvw / (caw + cvw + cgw) * dqsath_dT
!
!     else !lai_sno + stem_sno is > 0._dp but less than minplant. Vegetation present but not affecting heat fluxes.Need a red_pet for GPP.
!             !assign it the ground value.
!
!             red_pet(pft) = max(0._dp,gpet)  !this translates into demand so should be positive.
!             underpet(pft) = 0._dp
!             Eg_v(pft) = 0._dp
!             can_evap(pft) = 0._dp
!             Hg_v(pft) = 0._dp
!             dEgdt_v(pft) = 0._dp
!             dHgdt_v(pft) = 0._dp
!             dLgdt_v(pft) = 0._dp
!
!     end if
! end do !pft loop

end subroutine latsensflux

!----------------------------------------------------------------------------------------------------------------

function esat(temp)

!Function to calculate saturation vapor pressure in water and ice
!From CLM formulation, table 5.2, after Flatau et al. 1992

use parametersmod, only : sp,dp,Tfreeze

implicit none

real(sp) :: esat  !saturation vapor pressure (Pa)
real(sp), intent(in) :: temp !temperature in K

real(sp), dimension(9) :: al !coefficients for liquid water
real(sp), dimension(9) :: ai !coefficients for ice

real(sp), dimension(0:8) :: a !coefficients

real(sp) :: T

integer :: i

al(1) = 6.11213476
al(2) = 4.44007856e-1
al(3) = 1.43064234e-2
al(4) = 2.64461437e-4
al(5) = 3.05903558e-6
al(6) = 1.96237241e-8
al(7) = 8.92344772e-11
al(8) =-3.73208410e-13
al(9) = 2.09339997e-16

ai(1) = 6.11123516
ai(2) = 5.03109514e-1
ai(3) = 1.88369801e-2
ai(4) = 4.20547422e-4
ai(5) = 6.14396778e-6
ai(6) = 6.02780717e-8
ai(7) = 3.87940929e-10
ai(8) = 1.49436277e-12
ai(9) = 2.62655803e-15

if (temp <= Tfreeze) then   !these coefficients are for temperature values in Celcius
  a(0:8) = ai
else
  a(0:8) = al
end if

T = temp - Tfreeze

if (abs(T) < 1e-2) T = 0.   ! Avoid underflow at exponential operation (Leo O Lai, Jun 2021)

esat = a(0)

do i = 1,8
  esat = esat + a(i) * T**i
end do

esat = 100. * esat

end function esat

!----------------------------------------------------------------------------------------------------------------

function desdT(temp)

!Function to calculate the first derivative of saturation vapor pressure in water and ice vs. temperature
!From CLM formulation, table 5.3, after Flatau et al. 1992

use parametersmod, only : sp,dp,Tfreeze

implicit none

real(sp) :: desdT    !derivative of saturation vapor pressure
real(sp), intent(in) :: temp !temperature in K

real(sp), dimension(9) :: bl !coefficients for liquid water
real(sp), dimension(9) :: bi !coefficients for ice

real(sp), dimension(0:8) :: b !coefficients

! real(sp) :: tmp

real(sp) :: T

integer :: i

bl(1) = 4.44017302e-1
bl(2) = 2.86064092e-2
bl(3) = 7.94683137e-4
bl(4) = 1.21211669e-5
bl(5) = 1.03354611e-7
bl(6) = 4.04125005e-10
bl(7) =-7.88037859e-13
bl(8) =-1.14596802e-14
bl(9) = 3.81294516e-17

bi(1) = 5.03277922e-1
bi(2) = 3.77289173e-2
bi(3) = 1.26801703e-3
bi(4) = 2.49468427e-5
bi(5) = 3.13703411e-7
bi(6) = 2.57180651e-9
bi(7) = 1.32268878e-11
bi(8) = 3.94116744e-14
bi(9) = 4.98070196e-17

if (temp <= Tfreeze) then
  b(0:8) = bi
else
  b(0:8) = bl
end if

T = temp - Tfreeze  !these coefficients are for temperature values in Celcius

if (abs(T) < 1e-2) T = 0.   ! Avoid underflow at exponential operation (Leo O Lai, Jun 2021)

desdT = b(0)

do i = 1,8
  desdT = desdT + b(i) * T**i
end do

desdT = 100. * desdT

end function desdT

end module latsensfluxmod
