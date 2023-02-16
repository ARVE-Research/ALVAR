module radfluxmod

! Module adapted from ARVE-DGVM to calculate soil physics (Leo Lai, Nov 2022)

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine radiativeflux(day,i,snl,dayhour,nighthour,dayl,lat,cldf,tdew,delta,sw,direct,diffuse)

! calculates the heat flux from the atmosphere to the land surface
! also the sunlight and shaded plant absorbed photosynthetically active
! radiation. Based upon CLM 4.0 and other sources as listed
! coded by JM May 9 08

! use arveparams,                 only : dp,Tfreeze,npft,band,sigsb,Esnow,Esoil,pi,minplant
! use statevars,                  only : sv,dt,dayl
! use soilstate_vars,             only : asurf_g,TairK,Swv,Swg_b,Swg_v,Lwv,Lwg_b,Lwg_v,dLgdT_b,&
!                                         Latm,Lwvd,Eveg,Egrnd,surf
! use metvarsmod,                 only : met,atm
! use pft_state_variables,        only : veg

use parametersmod, only : i4,sp,dp,pi,daysec,Tfreeze,d2r,r2d
use pftparmod,     only : npft
use statevarsmod,  only : band

implicit none

integer(i4)                  , intent(in) :: day
integer(i4)                  , intent(in) :: i
integer(i4)                  , intent(in) :: snl                  !index value of top snow layer (negative is more layers)
integer(i4)                  , intent(in) :: dayhour
integer(i4)                  , intent(in) :: nighthour
real(sp)                     , intent(in) :: dayl
real(dp)                     , intent(in) :: lat
real(sp)                     , intent(in) :: cldf                 !cloud cover fraction (1 is fully overcast)
real(sp)                     , intent(in) :: tdew                 !air dew point temperature
real(sp)                     , intent(in) :: delta
real(sp)                     , intent(in) :: sw                   !downwelling shortwave (W m-2)
real(sp)                     , intent(in) :: direct                   !downwelling shortwave (W m-2)
real(sp)                     , intent(in) :: diffuse                   !downwelling shortwave (W m-2)

! pointers
! logical, pointer                        :: polar                !true if polar night
! real(dp), pointer                       :: cldf                 !cloud cover fraction (1 is fully overcast)
! integer,  pointer                       :: snl                  !index value of top snow layer (negative is more layers)
! real(dp), pointer                       :: rZ0                  !zenith angle at solar noon (radians)
! real(dp), pointer                       :: zen                  !cosine of the solar zenith angle (rad)
! real(dp), pointer                       :: Tdew                 !air dew point temperature
! real(dp), pointer                       :: absrad_annual        !annual total radiation absorbed (J m-2)
! real(dp), pointer, dimension(:)         :: sw                   !downwelling shortwave (W m-2)
! real(dp), pointer, dimension(:,:)       :: temp_veg             !temperature of the vegetation (K)
! real(dp), pointer, dimension(:,:)       :: fabi                 !diffuse flux absorbed by vegetation
! real(dp), pointer, dimension(:,:)       :: fabd                 !direct flux absorbed by vegetation
! real(dp), pointer, dimension(:,:)       :: fabi_vir             !diffuse flux absorbed by vegetation (virtual leaves)
! real(dp), pointer, dimension(:,:)       :: fabd_vir             !direct flux absorbed by vegetation (virtual leaves)
! real(dp), pointer, dimension(:,:)       :: ftid                 !downward diffuse fluxes per unit incident direct beam
! real(dp), pointer, dimension(:,:)       :: ftii                 !downward diffuse fluxes per unit incident diffuse beam
! real(dp), pointer, dimension(:)         :: Tsoil                !soil temperature (K)
! real(dp), pointer, dimension(:)         :: Tsoiln               !soil temperature (K) previous time step
! real(dp), pointer, dimension(:,:)       :: ftdd                 !down direct flux below veg per unit dir flx
! real(dp), pointer, dimension(:)         :: lai_sno              !instantaneous leaf area index per pft(m2 m-2)
! real(dp), pointer, dimension(:)         :: lai_vir              !virtual leaf area index per pft(m2 m-2) (for phenology)
! real(dp), pointer, dimension(:)         :: lai_sun              !sunlight leaf area index per pft (m2 m-2)
! real(dp), pointer, dimension(:)         :: lai_sha              !shaded leaf area index per pft (m2 m-2)
! real(dp), pointer, dimension(:)         :: stem_sno             !instantaneous stem area index per pft(m2 m-2)
! real(dp), pointer, dimension(:)         :: APAR_sun             !visible only shortwave rad absorbed by vegetation averaged for leaf area (W m-2) sunlight
! real(dp), pointer, dimension(:)         :: APAR_sha             !visible only shortwave rad absorbed by vegetation averaged for leaf area (W m-2) shaded
! real(dp), pointer, dimension(:,:)       :: scatcf_tot           !fraction of intercepted radiation that is scattered (0 to 1)
! real(dp), pointer, dimension(:,:)       :: scatcf_tot_vir       !fraction of intercepted radiation that is scattered (0 to 1)(virtual leaves)
! real(dp), pointer, dimension(:)         :: fsun                 !sunlit fraction of canopy
! real(dp), pointer, dimension(:)         :: fsha                 !shaded fraction of canopy
! real(dp), pointer, dimension(:)         :: bigK                 !the optical depth of direct beam per unit leaf and stem area
! logical, pointer, dimension(:)          :: virtual_leaf         !true if the leaf being calculated is a virtual one to determine phenology.
! real(dp), pointer                       :: fsnow                !fraction of the gridcell covered by snow (fraction)

! local variables
real(sp) :: rlat
real(sp) :: rdelta
real(sp) :: sinlat
real(sp) :: sindel
real(sp) :: coslat
real(sp) :: cosdel
real(sp) :: t1
real(sp) :: tinv
real(sp) :: Zn
real(sp) :: Z0
real(sp) :: rZ0
real(sp) :: rZn

real(sp) :: dt
real(sp) :: dtime
real(sp) :: sunrise
real(sp) :: hdl
real(sp) :: zed
real(sp) :: zen
real(sp) :: calday
real(sp) :: phi
integer(i4) :: hr
integer(i4) :: h


integer(i4) :: k
integer(i4) :: pft                                !counter
real(sp), dimension(band,npft) :: Svt           !temporary variable
real(sp), dimension(band) :: Sgt_v,Swg          !SW absorbed by ground below vegetation, and bare ground
real(sp) :: plantLS                             !total LAI + stem indices.
real(sp), dimension(band) :: Swdir              !incident direct beam solar fluxes (W m-2)
real(sp), dimension(band) :: Swdif              !incident diffuse beam solar fluxes (W m-2)
real(sp) :: D                                   !dew point depression (Tdew-TairK)
real(sp) :: Sw_noondir                          !solar noon downwelling shortwave direct (W m-2)
real(sp) :: Sw_noondif                          !solar noon downwelling shortwave diffuse (W m-2)
real(sp) :: Sw_dir_i                            !downwelling shortwave direct (W m-2) for this short timestep
real(sp) :: Sw_dif_i                            !downwelling shortwave diffuse (W m-2) for this short timestep
real(sp) :: phi_dirbeam_dir                     !unscattered direct beam absorbed by the canopy (visible)
real(sp) :: phi_dirbeam_dif                     !scattered direct beam absorbed as diffuse radiation by the canopy (visible)
real(sp) :: phi_difbeam                         !incoming diffuse radiation absorbed by the canopy (visible)

real(sp) :: sumdir,sumdif

!parameters
real(dp), parameter :: mu_avg = 1._dp            !average inverse optical depth for longwave rad.
real(dp), parameter :: la = 10.77_dp             !coefficients for the correction of downwelling longwave dependent on cld cover
real(dp), parameter :: lb =  2.34_dp             !
real(dp), parameter :: lc =-18.44_dp             !
real(dp), parameter :: visnir_part = 0.5_dp      !partition of visible and near infrared (0.5 and 0.5)

!--------

! begin calculations
! Zenith angle calculations taken from airmass.f90 and diszen.f90 in ARVE-DGVM (Leo Lai, Nov 2022)
rlat   = d2r * lat
rdelta = d2r * delta

sinlat = sin(rlat)
sindel = sin(rdelta)
coslat = cos(rlat)
cosdel = cos(rdelta)

!------

!Eqn. 2.5
if (abs(lat - delta) < 90.0 .and. abs(lat + delta) >= 90.0) then
 t1 = 12.0
else
 t1 = (12.0 / pi) * acos(-tan(rlat) * tan(rdelta))
end if

tinv = 1.0 / t1

!Eqn. 2.9
if (abs(lat + delta) >= 90.0) then
  Zn = acos(sin(rlat) * sin(rdelta) - cos(rlat) * cos(rdelta)) / d2r
else
  Zn = 90.0
end if

!Eqn. 2.10
if (abs(lat - delta) >= 90.0) then
  Z0 = 90.0
else
  Z0 = abs(lat - delta)
end if

rZ0 = Z0 * d2r  !convert to radians
rZn = Zn * d2r

!-------------------------

if (i == 1) then
  dtime = dayl * 3600.              ! Day
else
  dtime = daysec - (dayl * 3600.)   ! Night
end if

if (i == 1) then            ! Day time
  hr = dayhour
else if (i == 2) then       ! Night time
  hr = nighthour
end if

dt = 3600.

!-------------------------

do h = 1, hr

  !find the what the fraction of day we are in this timestep.
  hdl = 0.5 * dtime / 3600.  !hrs
  sunrise = 12. - hdl
  ! zed = sunrise + (dt * (0.5 + real(counter - 1)) / 3600.)
  zed = sunrise + h + 0.5

  calday = real(day) + zed / 24._dp !calendar day plus the fraction of the present day.

  !Solar declination in radians:
  rdelta = delta * d2r  !NOTE since I changed this to the Berger delta, we lose the ability to have a changing delta as the day goes on (yet gain
                        !the ability to accurately do paeloruns). This is likely a very small change through the day.JM Oct 30 08

    !below is a different calculation for delta. However this calculation is only valid for modern simulations.
    !find the Earth orbit seasonal angle in radians
    !theta = 2._dp * pi * calday / dayspy
    !delta = 0.006918_dp - 0.399912_dp * cos(theta) + 0.070257_dp * sin(theta) - &
    !        0.006758_dp * cos(2._dp * theta) + 0.000907_dp * sin(2._dp * theta) - &
    !        0.002697_dp * cos(3._dp * theta) + 0.001480_dp * sin(3._dp * theta)

  ! Calday is the calender day for Greenwich, including fraction
  ! of day; the fraction of the day represents a local time at
  ! Greenwich; to adjust this to produce a true instantaneous time
  ! For other longitudes, we must correct for the local time change:
  ! local time based on the longitude and day of year
  ! then compute the local cosine solar zenith angle
  ! HOWEVER, we do not treat individual timezones discretely In ARVE, we keep only day/night timesteps
  ! thus day is day at any point on the globe (not moving as is the case in reality). We do not need to
  ! correct for the longitude as a result.

  phi = calday !+ (real(lon(i)-1)/real(plon))

  zen = sin(rlat) * sin(rdelta) - cos(rlat) * cos(rdelta) * cos(2._dp * pi * phi)

  !-------------------------


  ! Initial default values

  ! plantLS = 0._dp
  ! APAR_sun(1:npft) = 0._dp
  ! APAR_sha(1:npft) = 0._dp
  ! Swv(1:npft) = 0._dp
  ! Swg_b = 0._dp
  ! Swg_v(1:npft) = 0._dp
  ! Eveg(1:npft) = 0._dp
  ! Lwv(1:npft) = 0._dp
  ! Lwvd(1:npft) = 0._dp
  ! Lwg_v(1:npft) = 0._dp

  ! This routine has been adapted from the CLM routine in parts with other references as noted.
  if (i == 1 .and. dayl > 0._dp) then !or polar !only do SW during day.

    ! Partition the total day downwelling shortwave direct and diffuse components into the amount delivered in this short timestep
    ! From Wang et al. 2002 Ecol. Modelling 00 1-14. JM Oct 3 08. Regressions from Figure 2.

    ! Find the downwelling shortwave at solar noon for both diffuse and direct components
    Sw_noondir = (direct/hr) * (2.28_dp - 1.1_dp * rZ0 + 0.8_dp * rZ0**2 - 0.23_dp * rZ0**3)
    Sw_noondif = (diffuse/hr) * (1.73_dp - 0.81_dp * rZ0 + 0.58_dp * rZ0**2 - 0.16_dp * rZ0**3)

    Sw_dir_i = Sw_noondir * cos((acos(zen) - rZ0) / (0.5_dp * pi - rZ0) * 0.5_dp * pi) * zen / cos(rZ0)  !Eqn 5
    Sw_dif_i = Sw_noondif * cos((acos(zen) - rZ0) / (0.5_dp * pi - rZ0) * 0.5_dp * pi)  !Eqn 7

    ! Split solar flux into vis and near infrared (it is 50% in each)
    Swdir(1) = visnir_part * Sw_dir_i
    Swdif(1) = visnir_part * Sw_dif_i

    Swdir(2) = visnir_part * Sw_dir_i
    Swdif(2) = visnir_part * Sw_dif_i

    sumdir = sumdir + Sw_dir_i
    sumdif = sumdif + Sw_dif_i

    ! print *, day, i, hr, hdl, sunrise, zed, rZ0,acos(zen), dayl,sumdir,sumdif, Sw_noondir, sw, direct, diffuse

    !increment the annual total radiation total
    ! absrad_annual = absrad_annual + (Sw_dir_i + Sw_dif_i) * dt

    ! Shortwave radiative flux to the soil and canopy surface

    ! Solar radiation is conserved as the sum of the diffuse and direct shortwave equaling Sv + Sg + the reflected solar
    ! radiation

    ! do pft = 1,npft
    !
    !   plantLS = lai_sno(pft) + stem_sno(pft)
    !
    !   ! Find the total solar radiation absorbed by vegetation and ground. Virtual leaves will not be considered here.
    !   if (plantLS > minplant) then !if there is 'real' vegetation present.
    !
    !     ! Total solar radiation absorbed by the vegetation (both direct and diffuse)
    !     Svt(1:band,pft) = max(0._dp, (Swdir(1:band) * fabd(1:band,pft) + Swdif(1:band) * fabi(1:band,pft)))   !eq. 4.3
    !     Swv(pft) = sum(Svt(1:band,pft))
    !
    !     ! Total solar radiation absorbed by the ground (under vegetation)
    !     Sgt_v(1:band) = max(0._dp,(Swdir(1:band) * ftdd(1:band,pft) * (1. - asurf_g(1,1:band)) + (Swdir(1:band)&
    !                    *  ftid(1:band,pft) + Swdif(1:band) * ftii(1:band,pft)) * (1. - asurf_g(2,1:band)))) !eq. 4.4
    !     Swg_v(pft) = sum(Sgt_v(1:band))
    !
    !   end if
    !
    !   ! Find the total absorbed photosynthetically active radiation for both real and virtual leaves
    !
    !   ! Store the radiation absorbed in visible for photosynthesis routine for 'real' leaves
    !   if (lai_sno(pft) > 0._dp) then
    !
    !     k = 1  !1 is the visible band
    !     phi_dirbeam_dir = Swdir(k) * (1._dp - exp(-bigK(pft)*plantLS)) * (1._dp - scatcf_tot(k,pft))!CLM 4.0 Eqn 4.9
    !     phi_dirbeam_dif = max(0._dp, (Swdir(k) * fabd(k,pft) - phi_dirbeam_dir))        !CLM 4.0 Eqn 4.10
    !     phi_difbeam     = Swdif(k) * fabi(k,pft)                                        !CLM 4.0 Eqn 4.11
    !
    !     ! visible only shortwave rad absorbed by vegetation SUNLIGHT leaves averaged for the ground area(W m-2)
    !     if (lai_sun(pft) > 0._dp)   APAR_sun(pft) = (phi_dirbeam_dir + phi_dirbeam_dif * fsun(pft) &
    !                                                 + phi_difbeam * fsun(pft)) &
    !                                                 * (lai_sno(pft) / plantLS) / lai_sun(pft) !CLM 4.12
    !
    !     ! visible only shortwave rad absorbed by vegetation SHADED leaves averaged for the ground area(W m-2)
    !     if (lai_sha(pft) > 0._dp) APAR_sha(pft) = (phi_dirbeam_dif * fsha(pft) &
    !                                               + phi_difbeam * fsha(pft)) &
    !                                               * (lai_sno(pft) / plantLS)/ lai_sha(pft) !CLM 4.13
    !
    !   ! Find APAR_sun for virtual leaves (for phenology), all else set to 0
    !   else if (virtual_leaf(pft)) then
    !
    !     k = 1  !1 is the visible band
    !     phi_dirbeam_dir = Swdir(k) * (1._dp - exp(-bigK(pft) * (lai_vir(pft) + stem_sno(pft)))) * (1._dp - scatcf_tot_vir(k,pft))!CLM 4.0 Eqn 4.9
    !     phi_dirbeam_dif = max(0._dp, (Swdir(k) * fabd_vir(k,pft) - phi_dirbeam_dir))        !CLM 4.0 Eqn 4.10
    !     phi_difbeam     = Swdif(k) * fabi_vir(k,pft)                                        !CLM 4.0 Eqn 4.11
    !
    !     ! visible only shortwave rad absorbed by vegetation SUNLIGHT leaves averaged for the ground area(W m-2)
    !     if (lai_vir(pft) > 0._dp) APAR_sun(pft) = (phi_dirbeam_dir + phi_dirbeam_dif * fsun(pft) &
    !                                              + phi_difbeam * fsun(pft)) *&
    !                                                (lai_vir(pft)/ (lai_vir(pft) + stem_sno(pft)))&
    !                                                / lai_vir(pft)! CLM 4.0 Eqn 4.12
    !
    !   end if  !lai if loop
    !
    ! end do  !pft loop
    !
    ! ! total solar radiation absorbed by the ground (bare ground, no canopy above)
    ! Swg(1:band) = max(0._dp, (Swdir(1:band) * (1. - asurf_g(1,1:band)) + Swdif(1:band) * (1. - asurf_g(2,1:band))))  !eq. 4.5
    ! Swg_b = sum(Swg(1:band))

  end if

  ! !--------
  ! ! Longwave fluxes
  !
  ! ! Downwelling atmospheric longwave flux
  !
  ! ! the following equations for downwelling longwave come from
  ! ! Josey et. al 2003, JGR Atmospheres 108, doi:10.1029/2002JC001418  eqn. 14
  !
  ! D = Tdew - TairK
  ! Latm = sigsb * (TairK + la * cldf**2 + lb * cldf + lc + 0.84_dp * (D + 4.01_dp))**4
  !
  ! ! calc ground emissivity
  ! Egrnd = (Esnow * fsnow) + (Esoil * (1. - fsnow))
  !
  !
  ! ! the following calculations are required each iteration
  !
  ! ! Upwelling longwave fluxes
  ! do pft = 1,npft
  !
  !   plantLS = lai_sno(pft) + stem_sno(pft)
  !
  !   if (plantLS > minplant) then !if there is vegetation present
  !
  !     ! calc veg emissivity
  !     Eveg(pft) = 1._dp - exp(-(plantLS) / mu_avg)  !CLM eqn 4.25
  !
  !     !find the downward longwave below vegetation
  !     Lwvd(pft) = ((1._dp - Eveg(pft)) * Latm) + Eveg(pft) * sigsb * (temp_veg(1,pft)**4) &
  !                 + 4._dp * Eveg(pft) * sigsb * (temp_veg(1,pft)**3) * (temp_veg(2,pft) - temp_veg(1,pft))   !4.17 (CLM 4.0 Eqn 4.21)
  !
  !     !net longwave flux for vegetation
  !     Lwv(pft) = (2._dp - Eveg(pft) * (1._dp - Egrnd)) * Eveg(pft) * sigsb * (temp_veg(2,pft)**4) &
  !                     - Eveg(pft) * Egrnd * sigsb * (Tsoiln(snl+1)**4) &
  !                     - Eveg(pft) * (1._dp + (1._dp - Egrnd) * (1._dp - Eveg(pft))) * Latm   !eq. 4.19 (CLM 4.0 Eqn 4.23)
  !
  !     !find the NET longwave flux for ground (positive towards atmosphere) (CLM 4.0 Eqn 4.22)
  !     Lwg_v(pft) = Egrnd * sigsb * Tsoiln(snl+1)**4 -    &                                            !upwelling longwave
  !                  Egrnd * Lwvd(pft) -  &                                                             !downwelling below veg
  !                  4._dp * Egrnd * sigsb * Tsoiln(snl+1)**3 * (Tsoil(snl+1) - Tsoiln(snl+1))         !rate of change
  !
  !   end if
  !
  ! end do
  !
  ! !Longwave fluxes on bare ground
  !  Lwg_b = Egrnd * sigsb * Tsoiln(snl+1)**4 -    &                                    !upwelling longwave
  !          Egrnd * Latm + &                                                                 !downwelling
  !          4._dp * Egrnd * sigsb * Tsoiln(snl+1)**3 * (Tsoil(snl+1) - Tsoiln(snl+1))         !rate of change
  !
  !  dLgdT_b = -4._dp * Egrnd * sigsb * Tsoiln(snl+1)**3
  !
  ! ! the net longwave flux for ground (both vegetated and not) is calculated in surfaceheatflux

end do




end subroutine radiativeflux





! subroutine surface_rad(i,lat,dayl,Ratm,Patm,delta,Pjj,dsw_t,prec,cldf,dir,dif,lastpet,srad,tmp)
!
! use parametersmod, only : i4,sp,dp,daysec,Pstd,pi
! use radiationmod, only : airmasspars,airmass,c00,c80,initairmass
!
! implicit none
!
! !arguments
! integer(i4), intent(in) :: i
! real(dp), intent(inout) :: lat        !atmospheric pressure (Pa)
! real(sp), intent(inout) :: dayl       !atmospheric pressure (Pa)
! real(sp), intent(inout) :: Ratm       !atmospheric pressure (Pa)
! real(sp), intent(inout) :: Patm       !atmospheric pressure (Pa)
! real(sp), intent(inout) :: delta       !atmospheric pressure (Pa)
! real(sp), intent(inout) :: Pjj        !precipitation equitability index
! real(sp), intent(inout) :: dsw_t      !top-of-atmosphere insolation (kJ m-2 d-1)
! real(sp), intent(inout) :: prec       !precipitation mm/day
! real(sp), intent(inout) :: cldf       !24 hour mean cloud cover fraction
! real(sp), intent(inout) :: dir
! real(sp), intent(inout) :: dif
! real(sp), intent(inout) :: lastpet    !potential evapotranspiration from last day (mm/timestep)
! real(sp), intent(inout) :: srad       !24 hour total downwelling shortwave (KJ m-2)
! real(sp), dimension(:), intent(inout) :: tmp  !monthly mean temp (C)
!
! ! real(dp), pointer :: Patm       !atmospheric pressure (Pa)
! ! real(dp), pointer :: Pjj        !precipitation equitability index
! ! real(dp), pointer :: mbar       !mean daily air mass
! ! real(dp), pointer :: mc         !airmass at medium cosine Z angle
! ! real(dp), pointer :: ml         !airmass at bottom-quarter cos Z angle
! ! real(dp), pointer :: mo         !airmass at max cosine zenith angle
! ! real(dp), pointer :: dsw_t      !top-of-atmosphere insolation (kJ m-2 d-1)
! ! real(dp), pointer :: prec       !precipitation mm/day
! ! real(dp), pointer :: cldf       !24 hour mean cloud cover fraction
! ! real(dp), pointer :: direct     !direct-beam downwelling shortwave (kJ m-2 d-1)
! ! real(dp), pointer :: diffuse    !diffuse downwelling shortwave (kJ m-2 d-1)
! ! real(dp), pointer :: lastpet    !potential evapotranspiration from last day (mm/timestep)
! ! logical, pointer :: polar       !polar night if true.
! ! real(dp), pointer :: dswb       !24 hour total downwelling shortwave (KJ m-2)
! ! real(dp), pointer, dimension(:) :: tmp  !monthly mean temp (C)
!
! !variables
! real(sp) :: dtime
! real(sp) :: mbar       !mean daily air mass
! real(sp) :: mc         !airmass at medium cosine Z angle
! real(sp) :: ml         !airmass at bottom-quarter cos Z angle
! real(sp) :: mo         !airmass at max cosine zenith angle
! real(sp) :: tau   !direct insolation atmospheric turbidity factor
! real(sp) :: zeta0 !diffuse insolation atmospheric turbidity factor
! real(sp) :: x     !tropics indicator (tropical = 1, north (lowest T below 10 C) = 0, otherwise btwn 0 & 1)
! real(sp) :: fm    !atmospheric transmittance function
! real(sp) :: p     !relative atmospheric pressure (1=sea level)
! real(sp) :: tcm   !minimum temperature in the year (used as tropics indicator)
! real(sp) :: sun   !bright sunshine duration fraction, n/N (fraction)
!
! real(sp) :: direct     !direct-beam downwelling shortwave (kJ m-2 d-1)
! real(sp) :: diffuse    !diffuse downwelling shortwave (kJ m-2 d-1)
! real(sp) :: dswb
!
! type(airmasspars) :: air
!
! !parameters
! real(sp), parameter :: kp  = 0.500_dp !links absorption coeff. to trans. coeff.
! real(sp), parameter :: kag = 3.300_dp
! real(sp), parameter :: kan = 2.320_dp
! real(sp), parameter :: kn  = 0.686_dp !cloud parameter
! real(sp), parameter :: ag = 0.17_dp   !Surface shortwave albedo (average=0.17)
!
! !---------------
!
! if (i == 1) then
!   dtime = dayl * 3600.              ! Day
! else
!   dtime = daysec - (dayl * 3600.)   ! Night
! end if
!
! call airmass(lat,delta,dayl,Ratm,air)
!
! mbar = air%mbar
! mc = air%mc
! ml = air%ml
! mo = air%mo
!
! if (i == 1) then  !daytime only calculation
!
!   p = min(1._dp,Patm / Pstd)
!
!   tcm = minval(tmp) !find min temp in the year (used as a tropics indicator)
!
!   if (tcm < 10.0_dp) then
!     x = 0.0_dp
!   else if (tcm > 20.0_dp) then !tropics
!     x = 1.0_dp
!   else
!     x = sin(pi / 2.0_dp * (tcm / 10.0_dp - 1.0_dp))
!   end if
!
!   !--
!   !sun is simply 1 - cloud cover
!   sun = 1._dp - cldf
!
!   !--
!   !find direct insolation atmospheric turbidity factor
!
!   tau = exp(-0.115_dp * p * ((2.15_dp - 0.713_dp * x + exp(-6.74_dp / (prec + 1._dp))) &
!                     * exp(0.0971_dp * lastpet) - 0.65_dp * (1._dp - x) * Pjj))  !Eqn. 4.1
!
!   !find atmospheric transmittance function
!   fm = 0.01452_dp * (mbar + ml) * exp(1.403_dp * tau) - 0.1528_dp * mo + mc + 0.48700_dp * (mc - ml) + 0.2323_dp   !Eqn. 2.4 2nd term
!
!   !direct beam downwelling shortwave (kJ m-2 d-1)
!   direct = sun * tau**kp * dsw_t * tau**fm   !Eqn. 2.4
!
!   !find diffuse insolation atmospheric turbidity factor
!   zeta0 = 0.503_dp * exp(-1.2_dp * p * exp(-0.633_dp / (prec + 1._dp) - 0.226_dp * lastpet)) &
!            * 3.3_dp**ag * 2.32_dp**(1._dp - sun) * (1._dp - 0.686_dp * (1._dp - sun))  !Eqn. 4.2
!
!   !diffuse downwelling shortwave (kJ m-2 d-1)
!   diffuse = zeta0 * kag**ag * kan**(1.0_dp - sun) * (1 - kn * (1.0_dp - sun)) * (tau**kp * dsw_t - direct)   !Eqn. 2.5
!
!   !for gridded data, assign the total radiation (in kJ m-2 d-1) to met (used by phenology)
!       dswb = diffuse + direct  !FLAG maybe not anymore, possibly extraneous JM 14.12.2010
!
!   !convert to W m-2 per timestep
!   ! if (polar) then
!   !         direct = direct * 1000. / daysec     !polar day spreads radiation over whole 24 hr period
!   !         diffuse = diffuse * 1000. / daysec
!   ! else !normal day
!   !         direct = direct * 1000. / dtime      !convert kJ m-2 day to W m-2 timestep-1
!   !         diffuse = diffuse * 1000. / dtime
!   ! end if
!
!   direct = direct * 1000. / dtime      !convert kJ m-2 day to W m-2 timestep-1
!   diffuse = diffuse * 1000. / dtime
!
!   !reset lastpet to 0
!   lastpet = 0._dp
!
! else  !night time set radiation to 0
!
!         diffuse = 0._dp                       !no radiation when night
!         direct = 0._dp
!         dswb = 0._dp
!
! end if
!
! print *, i, dayl, dswb* 1000. / dtime , srad* 1000. / dtime , direct, dir* 1000. / dtime  , diffuse, dif* 1000. / dtime
!
! end subroutine surface_rad




































end module radfluxmod
