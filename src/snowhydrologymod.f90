module snowhydrologymod

implicit none

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine newsnowflux(i,dayl,present,snl,Wsno,zsno,tausno,temp,dprec,dz,zpos,zipos,Wliq,Wice,Tsoil,Tsoiln, &
                       fpc_grid,xs_can_drip,fsnow,qliq,qsno,qgrnd_l,qgrnd_s)

! Subroutine to calculate current day new snow flux, called before soil temperature routine

use parametersmod, only : i4,sp,dp,Tfreeze,daysec
use pftparmod,     only : npft
use statevarsmod,  only : nl,ns
! use canopywater, only : canopy_interception

implicit none

integer(i4),                  intent(in)    :: i
real(sp),                     intent(in)    :: dayl                ! Daylength (h)
logical,  dimension(:),       intent(in)    :: present             ! True if pft is present in gridcell
integer(i4),                  intent(inout) :: snl                 ! Index value of top snow layer (negative is more layers)
real(sp),                     intent(inout) :: Wsno                ! Snow water equivalent of the snowpack (mm)
real(sp),                     intent(inout) :: zsno                ! Total thickness of the snowpack (m)
real(sp),                     intent(inout) :: tausno              ! Non-dimensional snow age
real(sp),                     intent(inout) :: temp                ! Air temperature (c)
real(sp),                     intent(inout) :: dprec               ! Total daily precipitation (mm)
real(sp), dimension(ns:nl),   intent(inout) :: dz                  ! Thickness of the soil layers (m)
real(sp), dimension(ns:nl),   intent(inout) :: zpos                ! z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(sp), dimension(ns:nl+1), intent(inout) :: zipos               ! z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(sp), dimension(ns:nl),   intent(inout) :: Wliq                ! Soil liquid water content at layer midpoint (mm)
real(sp), dimension(ns:nl),   intent(inout) :: Wice                ! Soil ice content at layer midpoint (mm)
real(sp), dimension(ns:nl),   intent(inout) :: Tsoil               ! Soil temperature (K)
real(sp), dimension(ns:nl),   intent(inout) :: Tsoiln              ! Soil temperature (K) previous time step
real(sp), dimension(:),       intent(inout) :: fpc_grid            ! Instantaneous cumulative fpc_grid per pft
real(sp), dimension(:),       intent(inout) :: xs_can_drip         ! Excess canopy water (mm)
real(sp),                     intent(inout) :: fsnow               ! Fraction of the gridcell covered by snow (fraction)
real(sp),                     intent(inout) :: qliq                ! Liquid precipitation (kg m-2 sec-1)
real(sp),                     intent(inout) :: qsno                ! Snowfall (kg m-2 sec -1)
real(sp), dimension(:),       intent(inout) :: qgrnd_l             ! Total rate of liquid precip reaching the ground under canopy(mm s-1)
real(sp), dimension(:),       intent(inout) :: qgrnd_s             ! Total rate of solid precip reaching the ground under canopy(mm s-1)

! Parameters
real(sp),    parameter :: pdrysnow = 50.0       ! Density of dry snow falling at under -15 deg C (kg m-3)
real(sp),    parameter :: pwetsnow = 169.15775  ! Density of wet snow falling at over 2 deg C (kg m-3) (50. + 1.7 * 17.**1.5, CLM eqn 7.18)
real(sp),    parameter :: Tc       = 2.5        ! Critical threshold temperature separating rain from snow (deg C)
real(sp),    parameter :: zsno_max = 4.0        ! Maximum snow depth (m)
real(sp),    parameter :: z0mg     = 0.01       ! Momentum roughness length for bare soil (m, CLM eqn 3.49)
real(sp),    parameter :: psno_new = 100.0      ! Density of new snow (kg m-3)
integer(i4), parameter :: m        = 1          ! Exponent for fsnow calc. Suggested to be 1 for global applications

!local variables
real(sp)    :: dtime           ! Current timestep (second)
real(sp)    :: TairK           ! Air temperature (K)
real(sp)    :: fpc_sum         ! Instnataneous sum of all PFT fpc in gridcell
real(sp)    :: Prate           ! Top of the canopy precipitation rate (kg m-2 sec-1)
real(sp)    :: Fpliq           ! Rain to snow transfer scalar (fraction)
real(sp)    :: pnsno = 0.      ! Density of newly fallen snow (kg m-3)
real(sp)    :: dzsno = 0.      ! Change in snow depth (m)
integer(i4) :: pft

!-------------------------

fsnow = 0.
TairK = temp + Tfreeze

if (i == 1) then
  dtime = dayl * 3600.              ! Day
else
  dtime = daysec - (dayl * 3600.)   ! Night
end if

! Find top of canopy precipitation rate
Prate = dprec / daysec    ! mm s-1

! Determine state of precipitation
if (TairK <= Tfreeze) then                ! Then the rain to snow transfer scalar is set to 0
  Fpliq = 0.0
else if (TairK <= Tfreeze + 2.) then      ! If it is above freezing but below the Tc then some melting
  Fpliq = max(0., -54.632 + 0.2 * TairK)
else if (TairK <= Tfreeze + Tc) then      ! Tc is critical temperature threshold seperating rain from snow,Tc = 2.5
  Fpliq = 0.4
else                                      ! The temperature is too high for snow so the Fpliq is set to 1
  Fpliq = 1.0
end if

! Partition precip into snow or rain
qliq = Prate * Fpliq
qsno = Prate * (1.0 - Fpliq)

!-------------------------
! Calculate canopy water (interception, canopy drip and throughfall) this will determine the amount of snow/rain to reach the ground
! do pft = 1,npft
!  if (present(pft)) then
!     call canopy_interception(j,pft)
!  else
!     qgrnd_l(pft) = 0.
!     qgrnd_s(pft) = 0.
!  end if
! end do
!
! ! Sum the throughfall and canopy drip across all pfts here
! qliq = (1. - fpc_sum) * qliq + sum((qgrnd_l(:) + xs_can_drip(:) * Fpliq / dtime)* fpc_grid(:))
! qsno = (1. - fpc_sum) * qsno + sum((qgrnd_s(:) + xs_can_drip(:) * Fpliq / dtime) * fpc_grid(:))

!-------------------------
! If snow exist, determine if the snow is wet or dry
if (qsno > 0.0) then

  ! Calculate the density of the NEW snow based upon air temperature
  if (TairK > Tfreeze + 2.0) then
    pnsno = pwetsnow
  else if (TairK <= Tfreeze + 2.0 .and. TairK > Tfreeze - 15.0) then
    pnsno = pdrysnow + 1.7 * (TairK - Tfreeze + 15.0) ** 1.5
  else
    pnsno = pdrysnow
  end if

  ! Update the change in snow depth, total thickness of the snowpack and water content
  dzsno = (qsno * dtime) / pnsno
  zsno  = zsno + dzsno

  ! Some areas have climate that does not allow summer melting of snow which results in the
  ! the build up of snow (i.e. glacier building). The snow can not presently be shifted away by wind
  ! so we just set a maximum snow depth. Any snow above that is assumed to be blown away (i.e. it is not
  ! added to soil). JM 21.04.2011
  if (zsno > zsno_max) then
     zsno = zsno_max
     qsno = 0.0
     dzsno = 0.0
  end if

  Wsno  = Wsno + qsno * dtime

  if (zsno >= 0.01) then

    if (snl+1 == 1) then

      ! Initialize a snow layer
      snl           = -1
      dz(snl+1)     = zsno
      zpos(snl+1)   = -0.5 * dz(snl+1)
      ! zipos(snl)    = -dz(snl+1)
      zipos(snl+1)    = -dz(snl+1)
      tausno        = 0.0
      Tsoil(snl+1)  = min(Tfreeze,TairK)
      Tsoiln(snl+1) = min(Tfreeze,TairK)  ! This is for radflux (which requires a previous timestep value)
      Wice(snl+1)   = Wsno
      Wliq(snl+1)   =  0.0

    else

      ! Update existing top snow layer
      Wice(snl+1) = Wice(snl+1) + qsno * dtime
      dz(snl+1)   = dz(snl+1) + dzsno
      zpos(snl+1) = zipos(snl+1) - 0.5 * dz(snl+1)
      ! zipos(snl)  = zipos(snl+1) - dz(snl+1)
      zipos(snl+1)  = zipos(snl+1) - dz(snl+1)

    end if !snow layer exists?

  end if !need to initialize a layer?

end if !new snow today?

!-------------------------
! Calculate fractional area of the gridcell covered by snow (Niu and Yang 2007)
if (Wsno > 0.) then
  fsnow = tanh(zsno / (2.5 * z0mg * ((min(Wsno/zsno, 800.)) / psno_new)**m))
end if

end subroutine newsnowflux

!---------------------------------------------------------------------

subroutine snowdynamics(i,dayl,snl,Wsno,Wsno_old,zsno,tausno,tausno_old,dz,zpos,zipos, &
                        Wliq,Wice,Tsoil,Tsoiln,fice0,ithaw,fsnow)

use parametersmod, only : i4,sp,dp,Tfreeze,daysec
use pftparmod,     only : npft
use statevarsmod,  only : nl,ns,snomax

implicit none

integer(i4),                  intent(in)    :: i
real(sp),                     intent(in)    :: dayl          ! Daylength (h)
integer(i4),                  intent(inout) :: snl           ! Index value of top snow layer (negative is more layers)
real(sp),                     intent(inout) :: Wsno          ! Snow water equivalent of the snowpack (mm)
real(sp),                     intent(inout) :: Wsno_old      ! Snow water equivalent of the snowpack (mm) from previous timestep
real(sp),                     intent(inout) :: zsno          ! Total thickness of the snowpack (m)
real(sp),                     intent(inout) :: tausno        ! Non-dimensional snow age
real(sp),                     intent(inout) :: tausno_old    ! Non-dimensional snow age of previous timestep
real(sp), dimension(ns:nl),   intent(inout) :: dz            ! Thickness of the soil layers (m)
real(sp), dimension(ns:nl),   intent(inout) :: zpos          ! z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(sp), dimension(ns:nl+1), intent(inout) :: zipos         ! z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(sp), dimension(ns:nl),   intent(inout) :: Wliq          ! Soil liquid water content at layer midpoint (mm)
real(sp), dimension(ns:nl),   intent(inout) :: Wice          ! Soil ice content at layer midpoint (mm)
real(sp), dimension(ns:nl),   intent(inout) :: Tsoil         ! Soil temperature (K)
real(sp), dimension(ns:nl),   intent(inout) :: Tsoiln        ! Soil temperature (K) from previous timestep
real(sp), dimension(ns:nl),   intent(inout) :: fice0         ! Layer ice fraction, previous timestep
integer(i4), dimension(ns:nl),   intent(inout) :: ithaw
real(sp),                     intent(inout) :: fsnow         ! Fraction of the gridcell covered by snow (fraction)

! Parameters
real(sp), parameter :: Cliq     = 4.18800e3   ! Heat capacity of liquid water  (J kg-1 K-1)
real(sp), parameter :: Cice     = 2.11727e3   ! Heat capacity of ice (typical) (J kg-1 K-1)
real(sp), parameter :: Lf       = 3.337e5     ! Water latent heat of fusion (J kg-1)
real(sp), parameter :: pliq     = 1000.0      ! Density of water (kg m-3)
real(sp), parameter :: pice     =  917.0      ! Density of ice (kg m-3)
real(sp), parameter :: pwetsnow = 169.15775   ! Density of wet snow falling at over 2 deg C (kg m-3) (50. + 1.7 * 17.**1.5, CLM eqn 7.18)
real(sp), parameter :: z0mg     = 0.01        ! Momentum roughness length for bare soil (m, CLM eqn 3.49)

!local variables
real(sp)    :: dtime           ! Current timestep (second)
integer(i4) :: l
logical :: combdivflag
real(sp), dimension(-snomax:0) :: Cr          ! Snow fractional compaction rate (s-1)
real(sp) :: Cr1                               ! Compaction rate due to destructive metamorphism (s-1)
real(sp) :: Cr2                               ! Compaction rate due to overburden (s-1)
real(sp) :: Cr3                               ! Compaction rate due to melting (s-1)
real(sp) :: c1                                ! Coefficients for calculating the compaction rates
real(sp) :: c2
real(sp), parameter :: c3 = 2.777e-6          ! (s-1)
real(sp), parameter :: c4 = 0.04           ! (K-1)
real(sp), parameter :: c5 = 0.08           ! (K-1)
real(sp), parameter :: c6 = 0.023          ! (m3 kg-1)
real(sp), parameter :: n0 = 9.e5              ! (kg s m-2)
real(sp), parameter :: tao_0 = 1.e-6          ! (s-1)
real(sp) :: Ps                                ! snow load pressure (?)
real(sp) :: Nu                                ! viscosity coefficient for overburden compaction (kg s m-2)
real(sp) :: r1                                ! for soil age calculations (grain growth
real(sp) :: r2                                ! effect of melt
real(sp) :: r3                                ! , soot and dirt)
real(sp) :: del_tausno                        ! increase in snow age this timestep
real(sp) :: snonew                            ! increase in snow amount since last timestep (mm)
integer(i4) :: nsno                           ! current number of snow layers
integer(i4) :: ll                             ! index of the lower snow layer used in division calculations
integer(i4) :: ul                             ! index of the upper snow layer used in division calculations
integer(i4) :: sl                             ! negative index value for snow layer division
integer(i4) :: ls                             ! index for layer shift down
real(sp), dimension(snomax) :: dzsno          ! snow layer thickness (m)
real(sp), dimension(snomax) :: dzsno_old      ! snow layer thickness from previous timestep (m)
real(sp), dimension(snomax) :: wlsno          ! Wliq for snow layers (m3 m-3)
real(sp), dimension(snomax) :: wisno          ! Wice for snow layers (m3 m-3)
real(sp), dimension(snomax) :: tsnow          ! snow layer temperature (K)
real(sp) :: xice
real(sp) :: xliq
real(sp) :: xfra
real(sp) :: xsno
real(sp) :: Wtot                              !total water mass of the snow layer (kg m-2)
real(sp) :: Tpor                              !air filled porosity of the snow layer (fraction)
real(sp), dimension(ns:nl) :: fice            !layer ice fraction, current timestep
real(sp), dimension(-snomax:0):: psno            !layer snow density (kg m-3)

!snow layer minimum and maximum thickness limits (m)
!note there is no limit to the thickness of the bottom layer
!layers get thinner towards the top, as for soil

real(sp), dimension(snomax),  parameter :: dzmin = [ 0.010, 0.015, 0.025, 0.055, 0.115 ] !min
real(sp), parameter :: psno_new        = 100.        ! density of new snow (kg m-3)
integer(i4), parameter :: m = 1                             ! exponent for fsnow calc. Suggested to be 1 for global applications

!--------------------------------------------

fsnow = 0.

if (i == 1) then
  dtime = dayl * 3600.              ! Day
else
  dtime = daysec - (dayl * 3600.)   ! Night
end if

! Begin calculations

! Snow layer compaction

        !write(*,*)'compact',dz(0),wice(0)

Ps = 0.
combdivflag = .false.

do l = snl+1, 0

  Wtot = Wice(l) + Wliq(l)
  Tpor = 1. - (Wice(l) / pice + Wliq(l) / pliq) / dz(l)

        !write(*,'(a,i4,5f12.4)')'start snow dyn',l,Tpor,(Wice(l) + Wliq(l)) / dz(l),Wice(l),Wliq(l),dz(l)

  !to account for the current layer add half of the ice and water contents of the layer being
  !compacted (CLM 4 Eqn 7.44 15.06.2010 JM)
  Ps = Ps + Wtot * 0.5

  !allow compaction only for non-saturated node and higher ice lens node.
  if (Tpor > 0.001 .and. Wice(l) > 0.1) then

    !compaction from metamorphism

    if (Wice(l) / dz(l) <= 100.) then  !kg m-3
      c1 = 1.
    else
      c1 = exp(-0.046 * (Wice(l) / dz(l) - 100.))
    end if

    if (Wliq(l) / dz(l) > 0.01) then
      c2 = 2.
    else
      c2 = 1.
    end if

    Cr1 = -c3 * c1 * c2 * exp(-c4 * (Tfreeze - Tsoil(l)))

    !compaction from overburden

    Nu = n0 * exp(c5 * (Tfreeze - Tsoil(l)) + c6 * (Wice(l) / dz(l)))

    Cr2 = -Ps / Nu

    !compaction from melting

    fice(l) = Wice(l) / (Wice(l) + Wliq(l))

    if (ithaw(l) == 1) then
      Cr3 = -1. / dtime * max(0., (fice0(l) - fice(l)) / fice0(l))
    else
      Cr3 = 0.
    end if

    !total compaction rate (s-1)

    Cr(l) = Cr1 + Cr2 + Cr3

    !update snow layer thickness

    dz(l) = dz(l) * (1. + Cr(l) * dtime)

                !write(*,'(a,i4,2f12.4)')'new layer thickness',l,dz(l),(Wice(l)+Wliq(l))/dz(l)

    !check if this makes the layer denser than ice
    psno(l) = (Wice(l) + Wliq(l)) / dz(l)

    if (psno(l) > pice) then !if so then set the thickness so that it is an ice layer.

        dz(l) = (Wice(l) + Wliq(l))/ pice

    end if

      !if the snow layer thickness is less than 0 then remove the layer.
      if (dz(l) < 0.) then
        dz(l) = 0.
        Wice(l) = 0.
        Wliq(l) = 0.
      end if

  else

  end if

  Ps = Ps + Wtot !sum(Wice(snl+1:l-1) + Wliq(snl+1:l-1))  !amount of ice and water in the layers above current layer (kg m-2)

end do

!-------------------------------
!snow layer combination

!first pass to determine if any layer is nearly melted
do l = snl+1,0

  if (Wice(l) <= 0.01) then

!write(*,'(a,i5,5f15.3)')'start compaction:layer combo ',l,Wice(l),Wliq(l),dz(l),(Wice(l) + Wliq(l)) / dz(l),Tsoil(l)-Tfreeze

    Wliq(l+1) = Wliq(l+1) + Wliq(l)
    Wice(l+1) = Wice(l+1) + Wice(l)

    Tsoil(l) = Tsoil(l+1)

    Wice(l)  = Wice(l+1)
    Wliq(l)  = Wliq(l+1)
    dz(l)    = dz(l+1)
    snl      = snl + 1

   if (dz(l) <= 0.) cycle

   !check if this makes either layer denser than ice
    psno(l) = (Wice(l) + Wliq(l)) / dz(l)

     if (psno(l) > pice) then !if so then set the thickness so that it is an ice layer.
      dz(l) = (Wice(l) + Wliq(l))/ pice
    end if

        !write(*,'(a,i5,5f15.3)')'end compaction:layer combo ',l,Wice(l),Wliq(l),dz(l),(Wice(l) + Wliq(l)) / dz(l),Tsoil(l)-Tfreeze

  end if
end do

if (snl == 0) then
  !no snow layers remain
  Wsno = 0.
  zsno = 0.
  return

else

  Wsno = sum(Wice(snl+1:0) + Wliq(snl+1:0))
  zsno = sum(dz(snl+1:0))

  if (zsno < 0.01) then       !total snowpack is too thin

    Wsno = sum(Wice(snl+1:0))
    Wliq(1) = Wliq(1) + sum(Wliq(snl+1:0))
    snl = 0
    return

  else

    !check that all remaining layers are at minimum thickness, or combine
    do l = snl+1,0

      sl = l - snl

               ! write(*,'(a,4i5,4f12.4)')'thickness',l,sl,snl,snl+1,dz(l),dzmin(sl),Tsoil(l)-Tfreeze,(Wice(l)+Wliq(l))/dz(l)

      if (dz(l) < dzmin(sl)) then      !layer too thin, combine it...
        if (l == snl+1) then          !surface layer, combine with layer below
          ul = l
          ll = l+1
        else if (l == 0) then        !bottom layer, combine with layer above
          ul = l-1
          ll = l
        else                          !combine with the thinnest neighboring layer
          if (dz(l+1) < dz(l-1)) then
            ul = l
            ll = l+1
          else
            ul = l-1
            ll = l
          end if
        end if

        call comb(dz(ll),Wliq(ll),Wice(ll),Tsoil(ll),dz(ul),Wliq(ul),Wice(ul),Tsoil(ul))
         combdivflag = .true.

        if (ll-1 > snl+1) then
          do ls = ll-1,snl+2,-1
            Tsoil(ls) = Tsoil(ls-1)
            Wice(ls)  = Wice(ls-1)
            Wliq(ls)  = Wliq(ls-1)
            dz(ls)    = dz(ls-1)
          end do
        end if

        snl = snl + 1
        if (snl >= -1) exit

      end if
    end do
  end if
end if

!---------------------------------------------------
!snow layer division

!There are four state variables that describe the snowpack:
!thickness (dz), temperature (Tsoil), water content (Wliq), ice content (Wice)

!The top snow layer is always the thinnest layer, thickening with depth.
!The snow layer against the ground (bottom layer) has index 0, with layers numbered negatively upwards.

!When performing layer division, we use temporary variables to avoid needing to renumber
!the layers each time a new one is formed.
!At the end of the routine, we translate all of the temporary variables back to the state variables.

!New temperature scheme from CLM 4.0 implimented. Without this temp scheme it was noticed that ARVE had
!temp spikes after a division. This new scheme will correct that problem. See CLM 4.0 pg. 134-135.

do l = snl+1,0

  sl = l - snl
  dzsno(sl) = dz(l)
  dzsno_old(sl) = dz(l)
  wlsno(sl) = Wliq(l)
  wisno(sl) = Wice(l)
  tsnow(sl) = Tsoil(l)

                !write(*,'(i5,a18,2i5,3f8.3)')-snl,'start division: layers, assign:',l,sl,dzsno(sl),Tsoil(l)-Tfreeze,(Wice(l) + Wliq(l)) / dz(l)

end do

nsno = abs(snl)

if (nsno == 1) then
  if (dzsno(1) > 0.03) then          !subdivide layer 2 from layer 1 (1,2)

    nsno = 2
    ul = 1
    ll = 2

    dzsno(1) = 0.5 * dzsno(1)
    wisno(1) = 0.5 * wisno(1)
    wlsno(1) = 0.5 * wlsno(1)

    dzsno(2) = dzsno(1)
    wisno(2) = wisno(1)
    wlsno(2) = wlsno(1)
    tsnow(2) = tsnow(1)

  end if

else if (nsno > 1) then

  if (dzsno(1) > 0.02) then

    xsno = dzsno(1) - 0.02  !drr

    xfra = xsno / dzsno(1)    !propor
    xice = xfra * wisno(1)    !zwice
    xliq = xfra * wlsno(1)    !zwliq

    xfra = 0.02 / dzsno(1)
    wisno(1) = xfra * wisno(1)
    wlsno(1) = xfra * wlsno(1)
    dzsno(1) = 0.02

    call comb(dzsno(2),wlsno(2),wisno(2),tsnow(2),xsno,xliq,xice,tsnow(1))
    combdivflag = .true.

    if (nsno <= 2 .and. dzsno(2) > 0.07) then !subdivide layer 3 from layer 2 (1,2,3)

      nsno = 3

      dzsno(2) = 0.5 * dzsno(2)
      wisno(2) = 0.5 * wisno(2)
      wlsno(2) = 0.5 * wlsno(2)

      dzsno(3) = dzsno(2)
      wisno(3) = wisno(2)
      wlsno(3) = wlsno(2)

      !adjust the new layer temperature as in CLM 4.0 Eqn. 7.59
      if (dzsno_old(1) + dzsno_old(2) > 0.) then
        tsnow(3) = tsnow(2) - ((tsnow(1) - tsnow(2)) / ((dzsno_old(1) + dzsno_old(2)) * 0.5)) * (dzsno(2) * 0.5)
      else
        tsnow(3) = tsnow(2)
      end if

      if (tsnow(3) >= Tfreeze) then  !Eqn 7.60
       tsnow(3) = tsnow(2)
      else
        if (dzsno_old(1) + dzsno_old(2) > 0.) then
          tsnow(2) = tsnow(2) - ((tsnow(1) - tsnow(2)) / ((dzsno_old(1) + dzsno_old(2)) * 0.5)) * (dzsno(2) * 0.5)
        end if
      end if

    end if
  end if
end if

if (nsno > 2) then
  if (dzsno(2) > 0.05) then

    xsno = dzsno(2) - 0.05
    xfra = xsno / dzsno(2)
    xice = xfra * wisno(2)
    xliq = xfra * wlsno(2)

    xfra = 0.05 / dzsno(2)
    wisno(2) = xfra * wisno(2)
    wlsno(2) = xfra * wlsno(2)
    dzsno(2) = 0.05

    call comb(dzsno(3),wlsno(3),wisno(3),tsnow(3),xsno,xliq,xice,tsnow(2))
    combdivflag = .true.

    if (nsno <= 3 .and. dzsno(3) > 0.18) then !subdivide layer 4 from layer 3 (1,2,3,4)

      nsno = 4

      dzsno(3) = 0.5 * dzsno(3)
      wisno(3) = 0.5 * wisno(3)
      wlsno(3) = 0.5 * wlsno(3)

      dzsno(4) = dzsno(3)
      wisno(4) = wisno(3)
      wlsno(4) = wlsno(3)

      !adjust the new layer temperature as in CLM 4.0 Eqn. 7.59
      if (dzsno_old(2) + dzsno_old(3) > 0.) then
        tsnow(4) = tsnow(3) - ((tsnow(2) - tsnow(3)) / ((dzsno_old(2) + dzsno_old(3)) * 0.5)) * (dzsno(3) * 0.5)
      else
        tsnow(4) = tsnow(3)
      end if

        if (tsnow(4) >= Tfreeze) then !Eqn 7.60
         tsnow(4) = tsnow(3)
        else

           if (dzsno_old(2) + dzsno_old(3) > 0.) then
             tsnow(3) = tsnow(3) - ((tsnow(2) - tsnow(3)) / ((dzsno_old(2) + dzsno_old(3)) * 0.5)) * (dzsno(3) * 0.5)
           end if

        end if

    end if
  end if
end if

if (nsno > 3) then
  if (dzsno(3) > 0.11) then

    xsno = dzsno(3) - 0.11
    xfra = xsno / dzsno(3)
    xice = xfra * wisno(3)
    xliq = xfra * wlsno(3)

    xfra = 0.11 / dzsno(3)
    wisno(3) = xfra * wisno(3)
    wlsno(3) = xfra * wlsno(3)
    dzsno(3) = 0.11

    call comb(dzsno(4),wlsno(4),wisno(4),tsnow(4),xsno,xliq,xice,tsnow(3))
    combdivflag = .true.

    if (nsno <= 4 .and. dzsno(4) > 0.41) then !subdivide layer 5 from layer 4 (1,2,3,4,5)

      nsno = 5

      dzsno(4) = 0.5 * dzsno(4)
      wisno(4) = 0.5 * wisno(4)
      wlsno(4) = 0.5 * wlsno(4)

      dzsno(5) = dzsno(4)
      wisno(5) = wisno(4)
      wlsno(5) = wlsno(4)

      !adjust the new layer temperature as in CLM 4.0 Eqn. 7.59 (modified for cases where the snow jumps to
      !five layers (result of using weather generator)
      if (dzsno_old(3) + dzsno_old(4) > 0.) then
        tsnow(5) = tsnow(4) - ((tsnow(3) - tsnow(4)) / ((dzsno_old(3) + dzsno_old(4)) * 0.5)) * (dzsno(4) * 0.5)
      else
        tsnow(5) = tsnow(4)
      end if

        if (tsnow(5) >= Tfreeze) then !Eqn 7.60
         tsnow(5) = tsnow(4)
        else

         if (dzsno_old(3) + dzsno_old(4) > 0.) then
           tsnow(4) = tsnow(4) - ((tsnow(3) - tsnow(4)) / ((dzsno_old(3) + dzsno_old(4)) * 0.5)) * (dzsno(4) * 0.5)
         end if

        end if

    end if
  end if
end if

if (nsno > 4) then
  if (dzsno(4) > 0.23) then

    xsno = dzsno(4) - 0.23
    xfra = xsno / dzsno(4)
    xice = xfra * wisno(4)
    xliq = xfra * wlsno(4)

    xfra = 0.23 / dzsno(4)
    wisno(4) = xfra * wisno(4)
    wlsno(4) = xfra * wlsno(4)
    dzsno(4) = 0.23

    call comb(dzsno(5),wlsno(5),wisno(5),tsnow(5),xsno,xliq,xice,tsnow(4))
    combdivflag = .true.

  end if
end if

!---------------------------------
!update state variables

        snl = -nsno

        do l = snl+1,0
          sl = l - snl
          dz(l)    = dzsno(sl)
          Wice(l)  = wisno(sl)
          Wliq(l)  = wlsno(sl)
          Tsoil(l) = tsnow(sl)

          if (Tsoiln(snl+1) == 0.) then
            Tsoiln(snl+1) = min(Tfreeze,Tsoil(snl+1))
          end if

          !check if any layer is denser than ice
          psno(l) = (Wice(l) + Wliq(l)) / dz(l)

          if (psno(l) > pice) then !if so then set the thickness so that it is an ice layer.

              dz(l) = (Wice(l) + Wliq(l))/ pice
          end if

        end do

        zsno = sum(dz(snl+1:0))

!recalculate layer depths and interfaces from the soil surface upwards
        do l = 0,snl+1,-1
          ! zpos(l) = zipos(l) - 0.5 * dz(l)
          ! zipos(l-1) = zipos(l) - dz(l)

          zpos(l) = zipos(l+1) - 0.5 * dz(l)
          zipos(l) = zipos(l+1) - dz(l)
        end do

!update snow age -added May 5 08 JM
        snonew = max(0.,(Wsno - Wsno_old))

        if (snonew > 10.) then !whole snow surface is renewed

                   tausno = 0.

            else !snow is aging

                !effect of grain growth due to vapour diffusion
                r1 = 5.e3 * (1. / Tfreeze - 1. / Tsoil(snl+1))

                !effect of near and at freezing of melt water
                r2 = min(0.,r1 * 10.)

                r1 = exp(r1)
                r2 = exp(r2)

                !effect of dirt and soot
                 r3 = 0.3

                !find the snow aging for this timestep
                if (Wsno > 0. .and. Wsno <= 800.) then
                  del_tausno = tao_0 * (r1 + r2 + r3) * dtime
                else
                  del_tausno = 0.
                end if

                tausno = max(0.,((tausno_old + del_tausno)*(1. - 0.1 * max(0.,(Wsno - Wsno_old)))))

        end if

        if (Wsno > 0.) then
          ! fractional area of the gridcell covered by snow (Niu and Yang 2007)
           fsnow = tanh(zsno / (2.5 * z0mg * ((min(Wsno/zsno, 800.)) / psno_new)**m))
        end if


        !if the snow layers are combined or subdivided and the number of snow layers is less than the max
           !set taosno to zero
        if (combdivflag .and. snl < snomax) tausno = 0.

        !set the previous timestep variables for next timestep
         tausno_old = tausno
          Wsno_old = Wsno


end subroutine snowdynamics

!--------------------------------------------------------------------------------------

subroutine comb(dz,Wliq,Wice,t,dz2,Wliq2,Wice2,t2)

use parametersmod, only : sp,Tfreeze

implicit none

!Combines characteristics of adjoining snow layers, operating on state vars: dz, t, wliq, wice.
!The combined temperature is based on the equation:
! the sum of the enthalpies of the two elements = that of the combined element.
!from CLM3

!arguments

real(sp), intent(inout) :: dz    ! nodal thickness of no. 1 element being combined [m]
real(sp), intent(inout) :: Wliq  ! liquid water of element 1
real(sp), intent(inout) :: Wice  ! ice of element 1 [kg/m2]
real(sp), intent(inout) :: t     ! nodal temperature of elment 1 [K]
real(sp), intent(in)    :: dz2   ! nodal thickness of no. 2 element being combined [m]
real(sp), intent(in)    :: Wliq2 ! liquid water of element 2 [kg/m2]
real(sp), intent(in)    :: Wice2 ! ice of element 2 [kg/m2]
real(sp), intent(in)    :: t2    ! nodal temperature of element 2 [K]

real(sp), parameter :: Cliq     = 4.18800e3   ! Heat capacity of liquid water  (J kg-1 K-1)
real(sp), parameter :: Cice     = 2.11727e3   ! Heat capacity of ice (typical) (J kg-1 K-1)
real(sp), parameter :: Lf       = 3.337e5     ! Water latent heat of fusion (J kg-1)
real(sp), parameter :: pice     =  917.0      ! Density of ice (kg m-3)

!local variables
real(sp) :: dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
real(sp) :: wliqc ! Combined liquid water [kg/m2]
real(sp) :: wicec ! Combined ice [kg/m2]
real(sp) :: Tc    ! Combined node temperature [K]
real(sp) :: h     ! enthalpy of element 1 [J/m2]
real(sp) :: h2    ! enthalpy of element 2 [J/m2]
real(sp) :: hc    ! temporary

!--------------

dzc = dz + dz2
Wicec = Wice + Wice2
Wliqc = Wliq + Wliq2

h  = (Cice * Wice  + Cliq * Wliq)  * (t  - Tfreeze) + Lf * Wliq
h2 = (Cice * Wice2 + Cliq * Wliq2) * (t2 - Tfreeze) + Lf * Wliq2

hc = h + h2

!This follows CLM 4.0 Eqn 7.55 to correct an error and ensure that enthalpy is always
!conserved during combination.
  Tc = Tfreeze + (hc - Lf * Wliqc) / (Cice * Wicec + Cliq * Wliqc)

dz   = dzc
Wice = Wicec
Wliq = Wliqc
t = Tc

end subroutine comb

















































end module snowhydrologymod
