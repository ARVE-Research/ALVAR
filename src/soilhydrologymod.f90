module soilhydrologymod

! Module adapted from LPJ-LMFire and ARVE-DGVM for soil hydrology

use parametersmod, only : i2,i4,sp,dp,missing_sp,Tfreeze

implicit none


! Parameters (copied from ARVE-DGVM, Leo Lai Aug 2021)
real(dp), parameter :: pcanmax = 0.1            ! max storage of water in canopy (mm)
real(dp), parameter :: Timp = 0.05_dp           ! volumetric liquid water content at which the soil is impermeable (fraction)
real(dp), parameter :: Pmax =-1.e8              ! maximum allowed value for soil matric potential (Pa)
real(dp), parameter :: psimax = -1.5e5          ! wilting point potential of leaves (CLM 8.11)
real(dp), parameter :: Esoil  = 0.96_dp         ! thermal emissivity of bare soil (fraction)
real(dp), parameter :: Ewater = 0.96_dp         ! thermal emissivity of water (fraction)
real(dp), parameter :: Esnow  = 0.97_dp         ! thermal emissivity of snow (fraction)
real(dp), parameter :: minsno = 0.05_dp         ! minimum snow depth to consider burial of vegetation (m)
real(dp), parameter :: Tstune = 0.34_dp         ! Tuning factor to turn first layer T into surface T (CLM parameterization, pg 88).
real(dp), parameter :: CNfac  = 0.5_dp          ! Crank Nicholson factor between 0 and 1
real(dp), parameter :: oneminusCNfac = 1._dp - CNfac
real(dp), parameter :: pdrysnow =  50.00000_dp  ! density of dry snow falling at under -15 deg C (kg m-3)
real(dp), parameter :: pwetsnow = 169.15775_dp  ! density of wet snow falling at over 2 deg C (kg m-3) (50._dp + 1.7_dp * 17._dp**1.5, CLM eqn 7.18)
real(dp), parameter :: Tc = 2.5_dp              ! critical threshold temperature separating rain from snow (deg C)
real(dp), parameter :: lb = 1.e-5               ! baseflow parameter (mm s-1) CLM eqn. 7.117
real(dp), parameter :: kd = 0.04_dp             ! saturated soil hydraulic cond. contributing to baseflow (mm s-1) CLM eqn. 7.118
real(dp), parameter :: Wpond = 10._dp           ! max. quantity of water for ponding (kg m-2)
real(dp), parameter :: alpha = 3._dp            ! adjustable scale dependent parameter (aquifer, soilwaterflux) CLM 4.0 p.155
real(sp), parameter :: ealpha = exp(-alpha)
real(dp), parameter :: fdecay = 2.5             ! decay factor (m-1)(aquifer and infiltration)
real(dp), parameter :: Smin  = 0.033_dp         ! irreducible water saturation of snow (fraction)
real(dp), parameter :: z0mg  = 0.01_dp          ! momentum roughness length for bare soil (m, CLM eqn 3.49)
real(dp), parameter :: wfact = 0.3_dp           ! weighted factor for topographic index (Eq. 7.53 CLM3.0, Leo Lai Aug 2021)
real(dp), parameter :: pliq  = 1000._dp         ! density of water (kg m-3)
real(dp), parameter :: pice  =  917._dp         !density of ice (kg m-3)
real(dp), parameter :: a = 0.053_dp             ! From Balland and Arp for dry soil thermal conductivity calc

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine soilwaterflux(i,snl,validcell,dayl,daet,sunrise,sunset,sunrise_n,dayhour,nighthour,hprec,hprec_n,dz,dzmm,zpos,zposmm,&
                        zipos,ziposmm,dayprec,nightprec,whc,Ksat,Tsat,Tpor,Ffrz,Psat,Bexp,aetsoil,aet,fpc_grid,rootfracl,&
                        Wliq,Wice,Tsoil,Tsoiln,Tliq,Tice,Psi,Psi_eq,Ku,zw,fsat,surf_runoff,surf_infl,Wliq_surf,dTliq_b)

! Adapted from ARVE-DGVM for ALVAR model by Leo O Lai (Aug, 2021)

use pftparmod,    only : npft
use statevarsmod, only : nl,ns
use utilitiesmod, only : tridiag

implicit none

integer(i4),                     intent(in)    :: i              ! i=1 (day) i-2 (night)
integer(i4),                     intent(in)    :: snl
logical,                         intent(in)    :: validcell
real(sp),                        intent(in)    :: dayl           ! Daylength (h)
real(sp),                        intent(in)    :: daet           ! Daily actual evapotranspiration (mm d-1)
integer(i4),                     intent(in)    :: sunrise        ! Current day sunrise hour
integer(i4),                     intent(in)    :: sunset         ! Current day sunset hour
integer(i4),                     intent(in)    :: sunrise_n      ! Next day sunrise hour
integer(i4),                     intent(in)    :: dayhour        ! Day time hours (h) --> from current day sunrise to sunset
integer(i4),                     intent(in)    :: nighthour      ! Night time hours (h) --> from current day sunset to next day sunrise
real(dp),    dimension(:),       intent(in)    :: hprec          ! Hourly precipitation (mm)
real(dp),    dimension(:),       intent(in)    :: hprec_n        ! Next day hourly precipitation (mm)
real(sp),    dimension(ns:nl),   intent(in)    :: dz             ! Snow/soil layer thickness (m)
real(sp),    dimension(ns:nl),   intent(in)    :: dzmm           ! Snow/soil layer thickness (mm)
real(sp),    dimension(ns:nl),   intent(in)    :: zpos           ! Midpoint z position (depth) of soil layer
real(sp),    dimension(ns:nl),   intent(in)    :: zposmm         ! zpos in mm
real(sp),    dimension(ns:nl+1), intent(in)    :: zipos          ! Snow/soil layer interface z position (depth), positive downwards (dim = nl+1, including surface and column bottom interface)
real(sp),    dimension(ns:nl+1), intent(in)    :: ziposmm        ! zipos in mm
real(sp),                        intent(inout) :: dayprec
real(sp),                        intent(inout) :: nightprec
real(sp),    dimension(:),       intent(inout) :: whc            ! Soil water holding capcity / available water content (mm)
real(sp),    dimension(ns:nl),       intent(inout) :: Ksat           ! Soil water saturated conductivity (mm s-1)
real(sp),    dimension(ns:nl),       intent(inout) :: Tsat           ! Soil water volumetric water content at saturation (fraction / m3 m-3)
real(sp),    dimension(ns:nl),       intent(inout) :: Tpor           ! Soil volumetric porosity (liquid minus ice content (or stones?)) (fraction)
real(sp),    dimension(nl),       intent(inout) :: Ffrz           ! Fractional impermeable area as a function of soil ice content at a layer (fraction)
real(sp),    dimension(ns:nl),       intent(inout) :: Psat           ! Soil water matric potential at saturation (mm)
real(sp),    dimension(ns:nl),       intent(inout) :: Bexp           ! Soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(sp),    dimension(:),       intent(in) :: aetsoil
real(sp),    dimension(:),       intent(in) :: aet
real(sp),    dimension(:),       intent(in) :: fpc_grid
real(sp),    dimension(:,:),     intent(in) :: rootfracl
real(sp),    dimension(ns:nl),       intent(inout) :: Wliq           ! Soil liquid water content at layer midpoint (mm)
real(sp),    dimension(ns:nl),       intent(inout) :: Wice           ! Soil ice content at layer midpoint (mm)
real(sp),    dimension(ns:nl),       intent(inout) :: Tsoil          ! Soil temperature (K)
real(sp),    dimension(ns:nl),       intent(inout) :: Tsoiln         ! Soil temperature for precious timestep (K)
real(sp),    dimension(ns:nl),       intent(inout) :: Tliq           ! Soil volumetric water content (fraction / m3 m-3) / Fractional soil water content
real(sp),    dimension(ns:nl),       intent(inout) :: Tice           ! Soil volumetric ice content (fraction / m3 m-3) / Fractional soil ice content
real(sp),    dimension(ns:nl),       intent(inout) :: Psi            ! Soil water potential
real(sp),    dimension(ns:nl),       intent(inout) :: Psi_eq         ! Restriction for min of soil potential (mm)
real(sp),    dimension(ns:nl),       intent(inout) :: Ku             ! Soil water instantaneous (unsaturated) conductivity across layer boundary (mm s-1)
real(sp),                        intent(inout) :: zw             ! Mean water table depth (m)
real(sp),                        intent(inout) :: fsat           ! Saturated fraction / partial contributing area (fraction)
real(sp),                        intent(inout) :: surf_runoff    ! Surface water runoff (mm)
real(sp),                        intent(inout) :: surf_infl      ! Surface water infiltration (mm)
real(sp),                        intent(inout) :: Wliq_surf
real(sp),                        intent(inout) :: dTliq_b

! Local variables
real(dp), allocatable, dimension(:) :: prec_tot     ! Current time-step total precipitation (day or night) (mm)
real(dp), allocatable, dimension(:) :: aet_tot      ! Current time-step total actual evaptranspiration (day or night) (mm)

real(dp) :: prec_tot0
real(dp) :: aet_tot0

real(sp) :: dtime                     ! Current day / night length (sec)
real(sp) :: dt                        ! Time-step length (seconds)

integer(i4) :: jint                   ! Index of the soil layer directly above the water table
real(sp)    :: zwmm                   ! Mean water table depth (mm)
real(sp)    :: dzmm_aq                ! Aqufier thickness

! Subroutine variables
real(sp) :: haet                      ! Hourly actual evapotranspiration (mm h-1)
real(sp) :: qover                     ! Surface liquid water runoff (mm s-1)
real(sp) :: qinfl                     ! Surface liquid water infiltration into top soil layer (mm s-1)
real(sp) :: qliq                      ! Liquid precipitation (variable from ARVE where snow is also considered)
real(sp) :: qliq0                     ! Liquid water reaching soil surface (mm s-1)
real(sp) :: qsdew                     ! Dew flux (mm s-1)
real(sp) :: qseva                     ! Soil (surface) evaporation flux (mm s-1)
real(sp) :: qsurf

real(sp) :: Seffpor                   ! Liquid water top soil layer relative to effective porosity adjusted for saturated frac
real(sp) :: varV                      ! Variable V from B5 in CLM 3.5 documentation
real(sp) :: qinflmax                  ! Maximum soil infiltration capacity (mm s-1)
real(sp) :: fsatfrac
real(sp) :: ffilled_pores
real(sp) :: funsat

real(sp), dimension(nl,npft) :: demandperlevel
real(sp), dimension(nl) :: rootfeff   ! Water extraction by plants per soil level (mm s-1)
real(sp), dimension(nl) :: vol_eq     ! Equilibrium volumetric water content (fraction)
real(sp)                :: tempi
real(sp)                :: temp0
real(sp)                :: voleq1

real(sp) :: s1
real(sp) :: s2
real(sp) :: s3

! Variables for tridiagonal system for water transfer between layers
real(sp) :: nterm                     ! Numerator term used in flux calculations
real(sp) :: dterm                     ! Denominator term used in flux calculations
real(sp) :: Fin                       ! Water flux into the soil layer
real(sp) :: Fout                      ! Water flux out of the soil layer
real(sp) :: ddFinTliq0                ! Derivative of the water flux into the layer with respect to theta up
real(sp) :: ddFinTliq1                ! Derivative of the water flux into the layer with respect to theta down
real(sp) :: ddFoutTliq1               ! Derivative of the water flux out of the layer with respect to theta up
real(sp) :: ddFoutTliq2               ! Derivative of the water flux out of the layer with respect to theta down
real(sp) :: dPsi_eq
real(sp) :: Psi1
real(sp) :: ddPsiTliq1
real(sp), dimension(nl) :: Sr         ! Soil wetness
real(sp), dimension(nl) :: ddKuTliq   ! Derivative of unsaturated conductivity with respect to theta
real(sp), dimension(nl) :: ddPsiTliq  ! dDrivative of matric potential with respect to theta

real(dp), dimension(nl+1) :: avect    ! Vectors for the tridiagonal solver (Refer to CLM Hydrology chapter)
real(dp), dimension(nl+1) :: bvect
real(dp), dimension(nl+1) :: cvect
real(dp), dimension(nl+1) :: rvect
real(dp), dimension(nl+1) :: dTliq    ! Change in volumetric soil liquid water over time-step dt (fraction)

real(sp), dimension(nl) :: T_wat      ! Total fractional soil water content (liquid + ice) (fraction)
real(sp), dimension(ns:0) :: snopor
real(sp) :: qliqsnoin
real(sp) :: qliqsnout


real(sp) :: fdec     ! Decay factor (m-1) for 0.5 degree grid, see CLM 4.0 pg. 135-136
real(sp) :: fmax     ! % of pixels in a grid cell whose topographic index is >= grid cell mean topographic index

integer(i4) :: l
integer(i4) :: dl
integer(i4) :: h
integer(i4) :: hr
integer(i4) :: pft

!-------------------------
! Compile met variables based on dayhour / nighthour

if (i == 1) then            ! Day time

  hr = dayhour
  dtime = dayl * 3600.

else if (i == 2) then       ! Night time

  hr = nighthour
  dtime = 86400. - (dayl * 3600.)

end if

allocate(prec_tot(hr))
allocate(aet_tot(hr))

! Calculate prec and aet total for current day/night timestep
if (i == 1) then            ! Day time

  prec_tot = hprec(sunrise:sunset-1)

  aet_tot = daet / real(dayhour)

  dayprec = sum(prec_tot)

  ! prec_tot = sum(hprec) * (dayhour / 24.) / dayhour

else if (i == 2) then       ! Night time

  prec_tot(1:(24-sunset+1))  = hprec(sunset:24)
  prec_tot((24-sunset+2):hr) = hprec_n(1:sunrise_n-1)

  aet_tot = 0.0

  nightprec = sum(prec_tot)

  ! prec_tot = sum(hprec) * (nighthour / 24.) / nighthour

end if

prec_tot0 = sum(prec_tot)
aet_tot0  = sum(aet_tot)

deallocate(prec_tot)
deallocate(aet_tot)

! prec_tot0 = sum(hprec)

qliq = prec_tot0 / dtime
! qliq = qliq * (1. - sum(fpc_grid,mask=tree))

!------------------------------------------
! Update the surface ice content, snow water content, and surface infiltration

! ziposmm(0:nl) = zipos(0:nl) * 1.e3        !interface depth in mm
! qseva_r = qseva_tot / dtime               !evaporative flux in mm/s
! qsdew_r = qsdew_tot / dtime               !dew flux in mm/s



! Water flux between snow layers
if (-snl > 0) then

  !-------------------------
  ! Update top layer water content for liquid precipitation, dew, frost and sublimation
  Wice(snl+1) = Wice(snl+1) !+ qfrost_tot - qsubl_tot

  !-------------------------
  if (Wice(snl+1) < 0._dp) then
    Wliq(snl+1) = Wliq(snl+1) + Wice(snl+1)
    Wice(snl+1) = 0._dp
  end if

  !-------------------------
  Wliq(snl+1) = max(0._dp, Wliq(snl+1) + qliq * dtime) !+ qsdew_tot - qseva_tot)

  !-------------------------
  ! Test whc of snow and transport water down snow column
  qliqsnoin = 0._dp

  do l = snl+1, 0
    Tice(l) = min(Wice(l) / (dz(l) * pice), 1._dp)
    snopor(l)  = 1._dp - Tice(l)
    Tliq(l) = min(Wliq(l) / (dz(l) * pliq), snopor(l))
  end do

  do l = snl+1, 0

    Wliq(l) = Wliq(l) + qliqsnoin

    if (l <= -1) then

      if (snopor(l) < Timp .or. snopor(l+1) < Timp) then                 !no flow
        qliqsnout = 0._dp
      else                                                               !water flows out of this layer
        qliqsnout = max(0._dp, (Tliq(l) - Smin * snopor(l)) * dz(l))
        qliqsnout = min(qliqsnout, (1._dp - Tice(l+1) - Tliq(l+1)) * dz(l+1))
      end if
    else
      qliqsnout = max(0._dp, (Tliq(l) - Smin * snopor(l)) * dz(l))
    end if

    qliqsnout = qliqsnout * pliq
    Wliq(l) = Wliq(l) - qliqsnout
    qliqsnoin = qliqsnout

  end do

  !-------------------------
  qliq = qliqsnout / dtime

  ! The total water flux to the soil cannot exceed the total water in the snowpack
  qliq = min(qliq,sum(Wliq(snl+1:0)+Wice(snl+1:0))/dtime)

end if

qliq0 = qliq   ! amount of liquid water reaching the soil surface

!------------------------------------------

! rootfeff = aet_tot0 / dtime
if (i == 1) then

  do pft = 1, npft
    demandperlevel(1:nl,pft) = aet(pft) * rootfracl(pft,1:nl)
  end do

  do l = 1, nl
    rootfeff(l) = sum(demandperlevel(l,1:npft)) / dtime
  end do

else

  rootfeff = 0.

end if

! Substitute surface evporation with actual evapotranspirtaiton for me as net surface water change
qseva = 0.        ! mm
qsdew = 0.        ! Assume to be zero for now


!-------------------------
! Time-step dt = 3600 for now as subroutine written for constant 1 hour time-step
! dt = 3600.

!------------------------------------------

! Calculate the number of iterations and timestep for water flux as
! the fraction of the infiltration rate over the porosity of the top soil layer
Tice(1) = min(Tsat(1), Wice(1) / (dz(1) * pice))
Tpor(1) = Tsat(1) - Tice(1)

! estimate the surface water flux as the net of rain and evaporation
! FLAG try adding dew as it was left out originally JM 26.06.2011
! qsurf = abs(qliq0 * dtime - qseva_tot + qsdew_tot)
qsurf = abs(qliq0 - qseva + qsdew) * dtime

! Calculate the wetness of the top soil layer to select an appropriate timestep for the hydrology calculations
Sr(1) = min(1._dp,(Tice(1) + Tliq(1)) / Tsat(1))

! The sub-daily timestep is given by the net surface flux over the top layer soil wetness
if (qsurf > 0._dp) then
 if (Sr(1) > 1.e-4) then
   dl = min(48,1 + nint(qsurf / Sr(1)))
 else
   dl = 48
 end if
else
 dl = 1
end if

dt = dtime / real(dl)  !sec

!-------------------------
! Run water model in hourly time-step based on length of dayhour / nighthour
hourloop : do h = 1, dl

  ! Saturated fraction (with fractional permeability added)
        ! NOTE fmax is the percent of pixels in a grid cell whose topographic index is larger than or equal to the grid cell
        ! mean topographic index. It is read-in in arve_dgvm.f90 into the iovariables soil structure.
        ! water table depth is needed here in meters.
  ! Eq. 7.54 in CLM3.0 (Leo Lai, Aug 2021)
  zw = zipos(nl+1)

  do l = 1, nl
    zw = zw - (Tliq(l) / Tsat(l)) * dz(l)
  end do
  if (zw < 0.) zw = 0

  ! zw = zipos(nl) - 1.0      ! Mean water table depth
  zwmm = zw * 1.e3

  ! Determine which soil layer corresponds to the one immediately above the present water table height.
  jint = nl
  do l = 2, nl
    if(zw <= zipos(l)) then
       jint = l-1
       exit
    end if
  end do

  ! Determine aqufier thickness
  if(jint < nl) then
    dzmm_aq = dzmm(nl)
  else
    dzmm_aq = zwmm - zposmm(nl)
  end if

  !------------------------
  ! Update surface water content
  Wliq_surf = max(0. , Wliq_surf + (qliq + qsdew - qseva) * dt)     !Eq. 7.25 (mm)

  Tice(1:nl) = min(Tsat(1:nl), Wice(1:nl) / (dz(1:nl) * pice))   ! Eq. 7.56
  Tpor(1:nl) = max(0.01, (Tsat(1:nl) - Tice(1:nl)))
  Tliq(1:nl) = min(Tpor(1:nl), Wliq(1:nl) / (dz(1:nl) * pliq))      ! Eq 7.57
  T_wat(1:nl) = Tice(1:nl) + Tliq(1:nl)

  !-------------------------

  ! Fractional permeability from Niu & Yang, Hydrometeorology 2006,v. 7,p. 937
  do l = 1, nl
    if (Wliq(l) > 0.) then
      Ffrz(l) = max(0.0, (exp(-alpha * (1.0 - 0.0)) - ealpha) / (1.0 - ealpha))
    else
      Ffrz(l) = max(0.0, (exp(-alpha * (1.0 - (Wice(l) / (Wice(l) + Wliq(l))))) - ealpha) / (1.0 - ealpha))      ! zero when noe ice / snow
    end if
  end do

  fmax = 0.9        ! not entire sure that this means... but it controls runoff (Leo Lai Jul 2021)
  fdec = 0.5        ! decay factor   ARVE-params

  fsat = (1.0 - Ffrz(1)) * fmax * exp(-0.5 * zw * fdec) + Ffrz(1)       ! Soil saturated fraction, one when no ice

  ! fsat = wfact * min(1.0, exp(-zw))       !Eq. 7.53 CLM3.0 (Leo Lai, Aug 2021)

  fsatfrac = fsat

  ! Maximum soil infiltration capacity
  ffilled_pores = max(0.01, (Tliq(1) / (max(Timp,(Tsat(1) - Tice(1))))))   !CLM 4.0 eq. 7.65
  funsat = max(0.01, (1.0 - fsatfrac))                                       !CLM 4.0 eq. 7.65
  Seffpor = max(0.0, ((ffilled_pores - fsatfrac) / funsat))                     !CLM 4.0 eq. 7.65

  varV = Bexp(1) * Psat(1) / (0.5 * dzmm(1))        !CLM 4.0 eq. 7.66 layer thickness in mm!

  qinflmax = Ksat(1) * (1.0 + varV * (Seffpor - 1.0)) !CLM 4.0 eq. 7.64

  !-------------------------

  ! Calculate surface runoff and infiltration
   if (Tpor(1) < Timp) then    !qover is liquid water surface runoff, qliq0 is liquid water reaching soil surface.
    qover = qliq0
   else
    qover = fsatfrac * qliq0 + (1.0 - fsatfrac) * max(0.0, qliq0 - qinflmax)
   end if

  ! Calculate how much water inflitrates the soil (assume no snow cover)
  qinfl = qliq0 - qover

  surf_runoff = qover * dtime
  surf_infl   = qinfl * dtime

  !-------------------------

  ! Calculate the equilibrium water content based on the water table depth (CLM 4 Eqns 7.120 - 7.122)
  do l = 1, nl

    if (zwmm < ziposmm(l)) then   !fully saturated when zw is less than the layer top

      vol_eq(l) = Tsat(l)

    ! Use the weighted average from the saturated part (depth > zw) and the equilibrium solution for the
    ! rest of the layer
    else if (zwmm < ziposmm(l+1) .and. zwmm > ziposmm(l)) then !CLM 4 Eqn 7.121

       ! Find the equilbrium volumetric water content for the unsaturated part of the layer (Eqn. 7.122)
      tempi = 1.0
      temp0 = ((Psat(l) - zwmm + ziposmm(l)) / Psat(l))**(1.0 - 1.0 / Bexp(l))
      voleq1 = Psat(l) * Tsat(l) / (1.0 - 1.0 / Bexp(l)) / (zwmm - ziposmm(l)) * (tempi - temp0)

      ! Find the equilbrium volumetric water content for the total layer (Eqn. 7.121)
      vol_eq(l) = (voleq1 * (zwmm - ziposmm(l)) + Tsat(l) * (ziposmm(l+1) - zwmm)) / (ziposmm(l+1) - ziposmm(l))
      vol_eq(l) = min(Tsat(l), vol_eq(l))
      vol_eq(l) = max(vol_eq(l), 0.0)

    else  ! Layers fully above the water table (zw) (CLM 4 Eqn 7.120)

      tempi = ((Psat(l) - zwmm + ziposmm(l+1)) / Psat(l))**(1.0 - 1.0 / Bexp(l))
      temp0 = ((Psat(l) - zwmm + ziposmm(l)) / Psat(l))**(1.0 - 1.0 / Bexp(l))
      vol_eq(l) = Psat(l) * Tsat(l) / (1.0 - 1.0 / Bexp(l)) / (ziposmm(l+1) - ziposmm(l)) * (tempi - temp0)
      vol_eq(l) = max(vol_eq(l), 0.0)
      vol_eq(l) = min(Tsat(l), vol_eq(l))

    end if

      Psi_eq(l) = Psat(l) * (max(vol_eq(l) / Tsat(l),0.01_dp))**(-Bexp(l))  !CLM 4 Eqn 7.125
      Psi_eq(l) = max(Pmax, Psi_eq(l))

  end do

  !-------------------------

  ! Hydraulic conductivity and soil matric potential and their derivatives
  ! Based upon Campbell 1974
  ! Similar to Eq. 7.70 in CLM3.0 (Leo Lai, Aug 2021)
  do l = 1, nl

    s1 = 0.5 * (T_wat(l) + T_wat(min(nl,l+1))) / (0.5 * (Tsat(l) + Tsat(min(nl,l+1))))
    s1 = min(1.0, s1)
    s2 = Ksat(l) * s1**(2.0 * Bexp(l) + 3.0)

    !CLM 4.0 Eqn. 7.80
    Ku(l) = (1.0 - 0.5 * (Ffrz(l) + Ffrz(min(nl,l+1)))) * s2  * s1

    !CLM 4.0 Eqn. 7.115
    ddKuTliq(l) = (1.0 - 0.5 * (Ffrz(l) + Ffrz(min(nl,l+1)))) * (2.0 * Bexp(l) + 3.0)&
                              * s2  * 0.5 / Tsat(l)

    ! Calc the soil wetness
    Sr(l) = min(1.0, ((Tliq(l) + Tice(l)) / Tsat(l)))
    Sr(l) = max(0.01, Sr(l))

    Psi(l) = Psat(l) * Sr(l)**(-Bexp(l))
    Psi(l) = max(Pmax, Psi(l))

    ddPsiTliq(l) = -Bexp(l) * Psi(l) / (Sr(l) * Tsat(l))

  end do

  !-------------------------------------------------
  ! Set up r, a, b, and c vectors for tridiagonal solution

  ! Node l=1 (top)
  l = 1

    Fin = qinfl
    dterm = zposmm(l+1) - zposmm(l)
    dPsi_eq = Psi_eq(l+1) - Psi_eq(l)
    nterm = (Psi(l+1) - Psi(l)) - dPsi_eq
    Fout = -Ku(l) * nterm / dterm
    ddFoutTliq1 = -(-Ku(l) * ddPsiTliq(l) + nterm * ddKuTliq(l)) / dterm
    ddFoutTliq2 = -( Ku(l) * ddPsiTliq(l+1) + nterm * ddKuTliq(l)) / dterm

    rvect(l) =  Fin - Fout - rootfeff(l)
    avect(l) =  0._dp
    bvect(l) =  dzmm(l) * (1._dp/dt) + ddFoutTliq1
    cvect(l) =  ddFoutTliq2

  !----------
  ! Middle soil layers (l = 2 to gnl-1)
  do l = 2, nl-1

    dterm = zposmm(l) - zposmm(l-1)
    dPsi_eq = Psi_eq(l) - Psi_eq(l-1)
    nterm = (Psi(l) - Psi(l-1)) - dPsi_eq
    Fin    = -Ku(l-1) * nterm / dterm
    ddFinTliq0 = -(-Ku(l-1) * ddPsiTliq(l-1) + nterm * ddKuTliq(l-1)) / dterm

    ddFinTliq1 = -( Ku(l-1) * ddPsiTliq(l)   + nterm * ddKuTliq(l-1)) / dterm
    dterm = zposmm(l+1) - zposmm(l)
    dPsi_eq = Psi_eq(l+1) - Psi_eq(l)
    nterm = (Psi(l+1) - Psi(l)) - dPsi_eq
    Fout = -Ku(l) * nterm / dterm
    ddFoutTliq1 = -(-Ku(l) * ddPsiTliq(l) + nterm * ddKuTliq(l)) / dterm
    ddFoutTliq2 = -( Ku(l) * ddPsiTliq(l+1) + nterm * ddKuTliq(l)) / dterm

    rvect(l) =  Fin - Fout - rootfeff(l)
    avect(l) = -ddFinTliq0
    bvect(l) =  dzmm(l) / dt - ddFinTliq1 + ddFoutTliq1
    cvect(l) =  ddFoutTliq2

  end do

  !----------
  ! Mode l=gnl (bottom)
  l = nl

  if (l > jint) then !water table is in soil column

    dterm = zposmm(l) - zposmm(l-1)
    dPsi_eq = Psi_eq(l) - Psi_eq(l-1)
    nterm = (Psi(l) - Psi(l-1)) - dPsi_eq
    dPsi_eq = Psi_eq(l) - Psi_eq(l-1)
    nterm = (Psi(l) - Psi(l-1)) - dPsi_eq
    Fin = -Ku(l-1) * nterm / dterm

    ddFinTliq0 = -(-Ku(l-1) * ddPsiTliq(l-1) + nterm * ddKuTliq(l-1)) / dterm
    ddFinTliq1 = -( Ku(l-1) * ddPsiTliq(l)   + nterm * ddKuTliq(l-1)) / dterm
    Fout =  0._dp
    ddFoutTliq1 = 0._dp

    rvect(l) =  Fin - Fout - rootfeff(l)
    avect(l) = -ddFinTliq0
    bvect(l) =  dzmm(l) / dt - ddFinTliq1 + ddFoutTliq1
    cvect(l) =  0._dp

    !Next set up aquifer layer; hydrologically inactive
    rvect(l+1) = 0._dp
    avect(l+1) = 0._dp
    bvect(l+1) = dzmm_aq / dt
    cvect(l+1) = 0._dp

  else ! water table is below soil column

    ! Compute aquifer soil moisture as average of layer gnl and saturation
    Sr(l) = max(0.5_dp * (1._dp + T_wat(l) / Tsat(l)), 0.01_dp)
    Sr(l) = min(1.0_dp, Sr(l))

    !Compute Psi for aquifer layer
    Psi1 = Psat(l) * Sr(l)**(-Bexp(l))
    Psi1 = max(Pmax, Psi1)

    !Compute ddPsiTliq for aquifer layer
    ddPsiTliq1 = -Bexp(l) * Psi1 / (Sr(l) * Tsat(l))

    !First set up bottom layer of soil column
    dterm = zposmm(l) - zposmm(l-1)
    dPsi_eq = Psi_eq(l) - Psi_eq(l-1)
    nterm = (Psi(l) - Psi(l-1)) - dPsi_eq
    Fin = -Ku(l-1) * nterm / dterm
    ddFinTliq0 = -(-Ku(l-1) * ddPsiTliq(l-1) + nterm * ddKuTliq(l-1)) / dterm
    ddFinTliq1 = -( Ku(l-1) * ddPsiTliq(l)   + nterm * ddKuTliq(l-1)) / dterm
    dterm = zposmm(l+1) - zposmm(l)
    dPsi_eq = Psi_eq(l+1) - Psi_eq(l)
    nterm = (Psi1 - Psi(l)) - dPsi_eq
    Fout = -Ku(l) * nterm / dterm
    ddFoutTliq1 = -(-Ku(l) * ddPsiTliq(l) + nterm * ddKuTliq(l)) / dterm
    ddFoutTliq2 = -( Ku(l) * ddPsiTliq1 + nterm * ddKuTliq(l)) / dterm

    rvect(l) =  Fin - Fout - rootfeff(l)
    avect(l) = -ddFinTliq0
    bvect(l) =  dzmm(l) / dt - ddFinTliq1 + ddFoutTliq1
    cvect(l) =  ddFoutTliq2

    !next set up aquifer layer; dterm/nterm unchanged, qin=qout
    Fin = Fout  !Fout from layer gnl
    ddFinTliq0 = -(-Ku(l) * ddPsiTliq(l) + nterm * ddKuTliq(l)) / dterm
    ddFinTliq1 = -( Ku(l) * ddPsiTliq1   + nterm * ddKuTliq(l)) / dterm
    Fout = 0._dp                  ! zero-flow bottom boundary condition
    ddFoutTliq1 = 0._dp          ! zero-flow bottom boundary condition

    rvect(l+1) =  Fin - Fout
    avect(l+1) = -ddFinTliq0
    bvect(l+1) =  dzmm_aq / dt - ddFinTliq1 + ddFoutTliq1
    cvect(l+1) =  0._dp

  end if

  ! Tridiagonal solver to calculate change in water content in each layer

  call tridiag(avect(1:nl+1),bvect(1:nl+1),cvect(1:nl+1),rvect(1:nl+1),dTliq(1:nl+1))

  !----------------------------

  ! Renew the mass of liquid water (add in change from this timestep)
  Wliq(1:nl) = Wliq(1:nl) + dTliq(1:nl) * dzmm(1:nl)

  Wliq = max(0., Wliq)

  Tliq(1:nl) = Wliq(1:nl) / (dz(1:nl) * pliq)

  call drainage(dt,dz,Ksat,Tsat,Bexp,Ku,zw,fsat,Wliq,Wice,Tliq,Tice)

  where (Wliq(1:nl) < 1.e-2)
    Wliq(1:nl) = 1.e-2
  end where

  where (Tliq(1:nl) < 1.e-3)
    Tliq(1:nl) = 1.e-3
  end where

end do hourloop

end subroutine soilwaterflux

!---------------------------------------------------------------------

subroutine drainage(dt,dz,Ksat,Tsat,Bexp,Ku,zw,fsat,Wliq,Wice,Tliq,Tice)

! Subroutine to calculate drainage from soil layer based on equations in CLM 3.0 (coded by Leo O Lai, Aug 2021)

use statevarsmod, only : nl,ns

implicit none

real(sp),                   intent(in)    :: dt
real(sp), dimension(ns:nl), intent(in)    :: dz
real(sp), dimension(ns:nl),     intent(in)    :: Ksat       ! Soil water saturated conductivity (mm s-1)
real(sp), dimension(ns:nl),     intent(in)    :: Tsat       ! Soil water volumetric water content at saturation (fraction / m3 m-3)
real(sp), dimension(ns:nl),     intent(in)    :: Bexp       ! Soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(sp), dimension(ns:nl),     intent(in)    :: Ku         ! Soil water instantaneous (unsaturated) conductivity across layer boundary (mm s-1)
real(sp),                   intent(in)    :: zw         ! Mean water table depth (m)
real(sp),                   intent(in)    :: fsat       ! Saturated fraction / partial contributing area
real(sp), dimension(ns:nl),     intent(inout) :: Wliq       ! Soil liquid water content at layer midpoint (mm)
real(sp), dimension(ns:nl),     intent(inout) :: Wice       ! Soil ice content at layer midpoint (mm)
real(sp), dimension(ns:nl),     intent(inout) :: Tliq       ! Soil volumetric water content (fraction / m3 m-3) / Fractional soil water content
real(sp), dimension(ns:nl),     intent(inout) :: Tice       ! Soil volumetric ice content (fraction / m3 m-3) / Fractional soil ice content

real(sp) :: qdrai
real(sp) :: qdrai_wet
real(sp) :: qdrai_dry
real(sp) :: wliq_exc
real(sp) :: wliq_def
real(sp) :: ddKuTliq_nl     ! Derivative of unsaturated conductivity with respect to theta for the lowest layer
real(sp) :: dTliq_nl
real(sp) :: wb
real(sp) :: wb_sum1
real(sp) :: wb_sum2
real(sp) :: Sr
real(sp) :: qdrai_bot

integer(i4) :: l

!------------------
! Remove excess water in each layer, i.e. when infiltration leads to oversaturation (Eq. 7.116 CLM3.0)

! do l = 1, nl
!
!   if (Tliq(l) / Tsat(l) > 1.0) then
!
!     wliq_exc = Wliq(l) - Tsat(l) * (dz(l) * pliq)
!
!     Wliq(l) = Wliq(l) - wliq_exc
!
!     Tliq(l) = Tsat(l)
!
!   end if
!
!   ! TESTING (Leo)
!   ! if (l < nl) Wliq(l+1) = Wliq(l+1) + wliq_exc
!   ! if (Tliq(l) / Tsat(l) > 1.0) then
!   !
!   !   wliq_exc = Wliq(l) - Tsat(l) * (dz(l) * pliq)
!   !
!   !   Wliq(l) = Wliq(l) - wliq_exc
!   !
!   !   Tliq(l) = Tsat(l)
!   !
!   ! end if
!
! end do

!------------------
! Calculate lateral drainage from saturation fraction (fsat) (Eq. 7.117 & 7.118 CLM3.0)

! wb_sum1 = 0.
! wb_sum2 = 0.
!
! do l = nl-3, nl-1
!
!   Sr = min(1.0, Tliq(l) / Tsat(l))
!
!   wb_sum1 = wb_sum1 + Sr * dzmm(l) * Ku(l)
!   wb_sum2 = wb_sum2 + dzmm(l) * Ku(l)
!
! end do
!
! wb = wb_sum1 / wb_sum2
!
! qdrai_dry = (1.0 - fsat) * kd * wb ** (2 * Bexp(1) + 3)
!
! qdrai_wet = fsat * lb * exp(-zw)
!
! !------------------
! ! Allocatble lateral drainage to layers depending on vertical diffuion (Eq. 7.123 CLM3.0)
! ! Oringal CLM3.0 has 10 layers and allocate lateral drainage to layers 6 to 9. Here we calculate lateral drainge for layer 3 to 5
! ! However, lateral drainage current not applicable because of no inter-gridcell interaction
!
! do l = nl-3, nl-1
!
!   Wliq(l) = Wliq(l) - dt * (qdrai_dry + qdrai_wet) * dzmm(l) * Ku(l) / wb_sum2
!
! end do

!------------------
! Calculate drainage from bottom of soil column (from the bottom layer) (Eq. 7.116 CLM3.0)

! Chnage in hydraulic conductivity due to change in liquid water content of bottom layer
ddKuTliq_nl = (2.*Bexp(nl) + 3.) * Ksat(nl) * (Tliq(nl) / Tsat(nl)) ** (2.*Bexp(nl) + 3.) &
              * (1 / Tsat(nl))

! Hydraulic conductivity (mm s-1) of bottom layer at current moisture level
dTliq_nl = (Ku(nl) * dt) / (dz(nl) * pliq)

Wliq(nl) = Wliq(nl) - (Ku(nl) + ddKuTliq_nl * dTliq_nl) * dt

!------------------

qdrai = qdrai_wet + qdrai_dry

qdrai = qdrai * dt

! Update fractional soil water content
Tliq(1:nl) = Wliq(1:nl) / (dz(1:nl) * pliq)

end subroutine drainage

!---------------------------------------------------------------------

! subroutine aquifier(grid,dt,qover)
!
! use statevarsmod,   only : soilvars
! use soilstatemod, only : nl
!
! integer(i4), intent(in) :: grid
! real(sp), intent(in) :: dt
! real(sp), intent(inout) :: qover
!
! ! Pointers
! real(sp), pointer  :: Waquif_a             !water in unconfined aquifer below soil (mm)
! real(sp), pointer  :: Waquif_t             !water in aquifer within soil column (mm).
! real(sp), pointer, dimension(:) :: Ksat    !soil water saturated conductivity at layer midpoint (mm s-1)
! real(sp), pointer, dimension(:) :: Psi     !soil water potential at layer midpoint (mm)
! real(sp), pointer, dimension(:) :: Psat    !soil water matric potential at saturation (mm)
! real(sp), pointer, dimension(:) :: Tice    !soil ice content (fraction)
! real(sp), pointer, dimension(:) :: Tsat    !soil water volumetric water content at saturation (fraction)
! real(sp), pointer, dimension(:) :: Tpor    !soil volumetric porosity (liquid minus ice content (or stones?)) (fraction)
! real(sp), pointer, dimension(:) :: Wliq    !soil liquid water content at layer midpoint (mm)
! real(sp), pointer, dimension(:) :: Wice    !soil ice content at layer midpoint (mm)
! real(sp), pointer, dimension(:) :: Bexp    !
! real(sp), pointer, dimension(:) :: Psi_eq  !
! real(sp), pointer, dimension(:) :: Ku      !soil water instantaneous (unsaturated) conductivity across layer boundary (mm s-1)
! ! real(sp), pointer               :: qover   !liquid water surface runoff (kg m-2 sec-1)
! real(sp), pointer               :: dTliq_b !dTliq value for level gnl+1
!
! logical :: peat
! integer :: jint
! real(sp) :: zw                   !mean water table depth (m)
! real(sp) :: zwmm
! real(sp), dimension(nl) :: zpos    !z coordinate of the middle of the soil layer relative to soil surface (m),positive downwards
! real(sp), dimension(nl) :: zposmm
! real(sp), dimension(nl) :: dz
! real(sp), dimension(nl) :: dzmm    !thickness of the soil layers (mm)
! real(sp), dimension(nl+1) :: zipos   !z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
! real(sp), dimension(nl+1) :: ziposmm
! real(sp) :: dzmm_aq !aquifer layer thickness (mm)
!
! !parameters
! real(sp), parameter :: Sy_min = 0.2                    !the fraction of water volume that can be drained by gravity in an unconfined aquifer (mineral soil)
!
! real(sp), parameter, dimension(3) :: Sy_peat = [ 0.655, 0.26, 0.125 ] !FLAG test 0.655, 0.26, 0.125
!
! real(sp), parameter :: aquif_max = 5000._sp        !max water in aquifer below soil profile (mm)
! real(sp), parameter :: qdraimax = 5.5e-3           !This new version is from CLM 4.0 tech note p.155
!                                                    !old=4.5e-4 !maximum subsurface drainage when the grid avg water table is 0 (kg m-1 s-1)
! real(sp), parameter :: Wlmin = 0.01_sp             !min water content per layer (mm)
!
! real(sp), parameter :: ealpha        = exp(-alpha)                          !for computational efficiency
! real(sp), parameter :: oneoverealpha = 1._sp / (1._sp - exp(-alpha))        !for computational efficiency
!
! !variables
! real(sp) :: qrecharge                      !recharge to the aquifer (mm s-1) positive when water enters aquifer
! real(sp) :: qdrai                          !liquid water draining from the soil (mm s-1)
! real(sp) :: fimp                           !fraction of impermeable area
! integer :: l,i                             !counter
! real(sp) :: Wexc                           !total soil column water surplus (mm s-1)
! real(sp) :: wh_zt                          !water head at the water table depth
! real(sp), dimension(nl) :: eff_porosity    !measure of how much pore space is available, adjusted for ice (fraction)
! real(sp) :: dzmmsum                        !soil thickness between water table and bottom of soil column (mm)
! real(sp) :: icefracsum                     !fraction of the soil between water table and bottom of soil column that contains ice
! real(sp) :: Sr                             !fractional soil moisture
! real(sp) :: s1                             !average fractional moisture between water table and layer Waquif_t
! real(sp) :: ka                             !unsaturated hydraulic conductivity
! real(sp) :: smp1                           !
! real(sp) :: wh                             !
! real(sp) :: ws                             !water used to fill soil air pores regardless of water content
! real(sp) :: Waquif_tsub                    !
! real(sp) :: xsi                            !
! real(sp) :: xs1,xs                         !
! real(sp) :: available_Wliq                 !
! real(sp) :: zw_old                         !previous timestep's water level
! real(sp) :: Sy
!
!
!
!
! Ksat => soilvars(grid)%Ksat
! Psi => soilvars(grid)%Psi
! Psat => soilvars(grid)%Psat
! Tice => soilvars(grid)%Tice
! Tsat => soilvars(grid)%Tsat
! Tpor => soilvars(grid)%Tpor
! Wliq  => soilvars(grid)%Wliq
! Wice  => soilvars(grid)%Wice
! Bexp => soilvars(grid)%Bexp
! Psi_eq => soilvars(grid)%Psi_eq
! Ku => soilvars(grid)%Ku
!
! Waquif_a => soilvars(grid)%Waquif_a
! Waquif_t => soilvars(grid)%Waquif_t
! dTliq_b => soilvars(grid)%dTliq_b
!
! peat = .false.
!
! ! Midpoint of soil layer zpos = [-2.5, -10, -22.5, -45, -80, -150] cm
! ! dz = soil layer depth in meter
! dz = [0.05, 0.1, 0.15, 0.3, 0.4, 1.0]           ! Soil layer thickness in meter
! dzmm = dz * 1.e3
!
! zpos = [0.025, 0.1, 0.225, 0.45, 0.8, 1.5]   ! midpoint of soil layer in meter, positive downwards
! zposmm = zpos * 1.e3
!
! zipos = [0., 0.05, 0.15, 0.3, 0.6, 1.0, 2.0]    ! soil position z position in meter, positive downwards (nl + 1, including surface and bottom interface)
! ziposmm = zipos * 1.e3
!
! ! Saturated fraction (with fractional permeability added)
!       ! NOTE fmax is the percent of pixels in a grid cell whose topographic index is larger than or equal to the grid cell
!       ! mean topographic index. It is read-in in arve_dgvm.f90 into the iovariables soil structure.
!       ! water table depth is needed here in meters.
! zw = zipos(nl) - 1.0      ! Mean water table depth
! zwmm = zw * 1.e3
!
! ! Determine which soil layer corresponds to the one immediately above the present water table height.
! jint = nl
! do l = 2, nl
!   if(zw <= zipos(l)) then
!      jint = l-1
!      exit
!   end if
! end do
!
! ! Determine aqufier thickness
! if(jint < nl) then
!   dzmm_aq = dzmm(nl)
! else
!   dzmm_aq = zwmm - zposmm(nl)
! end if
!
!
!
!
!
!
!
!
! !--------------------
!
! !Initial presets.
!   Wexc = 0._dp
!   zw_old = zw
!
!   if (.not. peat) then
!   Sy = Sy_min
!   else !peat
!       if (zpos(jint) < 0.3) then
!       Sy = Sy_peat(1)
!       else if (zpos(jint) > 0.3 .and. zpos(jint) < 1.) then
!       Sy = Sy_peat(2)
!       else
!       Sy = Sy_peat(3)
!       end if
!   end if
!
!        ! calculate qrecharge
!        if(jint < nl) then !if the water table is within the soil column
!
!           !water head at the water table depth
!           wh_zt = 0._dp   !defined as zero.
!
!           !Sr = max(Wliq(jint)/dzmm(jint)/Tsat(jint), 0.01_dp)
!           Sr = max(Wliq(jint)/dzmm(jint)/(Tsat(jint) - Wice(jint)/dzmm(jint)), 0.01_dp) !FLAG test
!           Sr = min(1._dp, Sr)
!
!           !use average moisture between water table and layer Waquif_t (which has to be 1)
!           s1 = 0.5_dp * (1._dp + Sr)
!           s1 = min(1._dp, s1)
!
!           !this is the expression for unsaturated hydraulic conductivity
!           ka = Ksat(jint) * s1**(2._dp * Bexp(jint) + 3._dp)
!
!           ! Recharge rate qrecharge to groundwater (positive to aquifer)
!           smp1 = Psat(jint) * Sr**(-Bexp(jint))
!           smp1 = max(Pmax,smp1)
!           wh   = smp1 - Psi_eq(jint)
!           qrecharge = -ka * (wh_zt - wh)  /((zw - zpos(jint)) * 1000._dp)
!
!        else
!
!           !if water table is below soil column, compute qrecharge from dTliq(gnl+1)
!           qrecharge = dTliq_b * dzmm_aq / dt
!
!        end if
!
!      ! To limit qrecharge  (for the first several timesteps)
!           qrecharge = max(-10.0_dp / dt,qrecharge)
!           qrecharge = min( 10.0_dp / dt,qrecharge)
!
!     eff_porosity(1:nl) = max(0.01_dp,Tsat(1:nl) - Tice(1:nl))
!
!     ! Topographic runoff
!
!        dzmmsum = sum(dzmm(jint:nl))
!
!        icefracsum = sum((Wice(jint:nl) / (Wice(jint:nl) + Wliq(jint:nl)) * dzmm(jint:nl)))
!
!        fimp = max(0._dp,(exp(-alpha * (1._dp - (icefracsum / dzmmsum))) - ealpha) * oneoverealpha)
!
!        qdrai = (1._dp - fimp) * qdraimax * exp(-fdecay * zw)
!
! if (peat) qdrai = 0.  !FLAG!! test
!
!     ! Water table calculation
!
!        ! Water storage in aquifer + soil
!        Waquif_t = Waquif_t + (qrecharge - qdrai) * dt
!
!        if (jint == nl) then             ! water table is below the soil column
!
!           Waquif_a  = Waquif_a + (qrecharge - qdrai) * dt
!           Waquif_t  = Waquif_a
!           zw = (zipos(nl) + 25._dp) - Waquif_a / 1000._dp / Sy
!           Wliq(nl) = Wliq(nl) + max(0._dp,(Waquif_a - aquif_max))
!           Waquif_a  = min(Waquif_a, aquif_max)
!
!        else   ! water table within soil layers
!
!           if (jint == nl-1) then       ! water table within bottom soil layer
!
!              zw = zipos(nl)- (Waquif_t - Sy * 1000._dp * 25._dp) / eff_porosity(nl) / 1000._dp
!
!           else  ! water table within soil layers 1-nl-1
!
!              ws = 0._dp   ! water used to fill soil air pores regardless of water content
!
!              do l = jint+2,nl
!                ws = ws + eff_porosity(l) * dzmm(l)
!              end do
!
!              zw = zipos(jint+1) - (Waquif_t - Sy * 1000_dp * 25._dp - ws) / eff_porosity(jint+1) / 1000._dp
!
!           end if
!
!           ! Remove subsurface runoff
!
! 	  ! NOTE! CLM subsurface runoff does not make sense. As the lateral subsruface runoff is removed
! 	  ! it simply disappears from the model as there is no lateral transfer of water between gridcells.
!           ! but excluding it leaves to instabilities. Necessary evil for the moment.JM 05.07.2011
!
!          if (.not. peat) then !FLAG!!
!
!              Waquif_tsub = 0._dp
!              do l = jint+1, nl
!                 Waquif_tsub = Waquif_tsub + Ku(l) * dzmm(l)
!              end do
!
!              do l = jint+1, nl
!                if (Waquif_tsub > 0._dp) then
!                    Wliq(l) = Wliq(l) - qdrai * dt * Ku(l) * dzmm(l) / Waquif_tsub
!                 end if
!              end do
!
!           end if
!
!        end if
!
!          !Check for large jumps in water table and limit
!         if (ABS(zw_old - zw) > 0.1) then
!
!           if (zw_old - zw > 0.d0) then
!               zw = zw_old - 0.1
!           else
!                 zw = zw_old + 0.1
!           end if
!
!         end if
!
!         zw = max(0.05_dp,zw)
!         zw = min(80._dp,zw)
!
!     ! Handle any excess water
!
!          !  excessive water above saturation added to the above unsaturated layer like a bucket
!          !  if column fully saturated, excess water goes to runoff
!
!          do l = nl,2,-1
!                xsi = max(Wliq(l) - eff_porosity(l) * dzmm(l),0._dp)
!                Wliq(l) = min(eff_porosity(l) * dzmm(l), Wliq(l))
!                Wliq(l-1) = Wliq(l-1) + xsi
!          end do
!
!          ! top layer
!             xs1 = max(max(Wliq(1),0._dp) - max(0._dp,(Wpond + Tsat(1) * dzmm(1) - Wice(1))),0._dp)
!             Wliq(1) = min(max(0._dp, Wpond + Tsat(1) * dzmm(1) - Wice(1)), Wliq(1))
!             Wexc     = xs1 / dt
!
!             ! add the excess water to the over land run-off
!             ! NOTE: this is a change from CLM which adds it to qdrai. However, that does not make
!             ! sense that the excess water at the surface should be added to subsurface drainage since it is not
!             ! below the surface. We add it to over land runoff instead. JM 27.10.10
!             qover = qover + Wexc
!
!     ! Handle any excess dryness
!
!            ! Limit Wliq to be greater than or equal to Wlmin.
!            ! Get water needed to bring Wliq equal Wlmin from lower layer.
!            ! If insufficient water in soil layers, get from aquifer water
!
!            do l = 1,nl-1
!                  if (Wliq(l) < Wlmin) then
!                     xs = Wlmin - Wliq(l)
!                  else
!                     xs = 0._dp
!                  end if
!                  Wliq(l) = Wliq(l) + xs
!                  Wliq(l+1) = Wliq(l+1) - xs
!            end do
!
!        ! Get water for bottom layer from layers above if possible
!            l = nl
!               if (Wliq(l) < Wlmin) then
!                  xs = Wlmin - Wliq(l)
!                  searchforwater: do i = nl-1,1,-1
!
!                     available_Wliq = max(Wliq(i) - Wlmin - xs,0._dp)
!
!                     if (available_Wliq >= xs) then
!                       Wliq(l) = Wliq(l) + xs
!                       Wliq(i) = Wliq(i) - xs
!                       xs = 0._dp
!                       exit searchforwater
!                     else
!                       Wliq(l) = Wliq(l) + available_Wliq
!                       Wliq(i) = Wliq(i) - available_Wliq
!                       xs = xs - available_Wliq
!                     end if
!                  end do searchforwater
!               else
!                  xs = 0._dp
!               end if
!
!        ! Needed in case there is no water to be found
!               Wliq(l) = Wliq(l) + xs
!               Waquif_t = Waquif_t - xs
!
!         dTliq_b = 0._dp
!
!
! end subroutine aquifier


end module soilhydrologymod
