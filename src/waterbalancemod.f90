module waterbalancemod

! Water balance module copied from LPJ-LMFire (Feb 2021 by Leo O Lai)
! Include alterantive simplesoilLPJ() subroutine for two layer soil model

implicit none

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine waterbalance(dpet,awc,swf,dgp,dphen,present,rm_ind,fpc_grid,rootfracl,daet,aet,aetsoil,dwscal,dgc)

use parametersmod, only : i4,sp
use pftparmod,     only : npft,pftpar
use statevarsmod,  only : nl,ns

implicit none

real(sp),                 intent(in)    :: dpet
real(sp), dimension(:),   intent(in)    :: awc
real(sp), dimension(nl),   intent(in)    :: swf
real(sp), dimension(:),   intent(in)    :: dgp      ! Non-water stressed canopy conductance (?)
real(sp), dimension(:),   intent(in)    :: dphen    ! Phenological state (0=no leaves, 1=full canopy)
logical,  dimension(:),   intent(in)    :: present
real(sp), dimension(:),   intent(in)    :: rm_ind
real(sp), dimension(:),   intent(in)    :: fpc_grid
real(sp), dimension(:,:), intent(in)    :: rootfracl
real(sp),                 intent(inout) :: daet
real(sp), dimension(:),   intent(inout) :: aet
real(sp), dimension(:),   intent(inout) :: aetsoil
real(sp), dimension(:),   intent(inout) :: dwscal
real(sp), dimension(:),   intent(inout) :: dgc      ! Actual canopy conductance (?)

! Parameters
real(sp), parameter :: emax   = 5.   ! Maximum daily transpiration rate (mm/day)
real(sp), parameter :: alpham = 1.4  ! Maximum Priestly-Taylor coefficient
real(sp), parameter :: gm     = 5.   ! Scaling conductance (mm s-1)?

! Local variables
real(sp), dimension(nl)   :: rootprop
real(sp), dimension(nl)   :: beta
real(sp) :: supply
real(sp) :: demand
real(sp) :: demandpot
real(sp) :: wr

integer(i4) :: pft

!-------------------------

daet     = 0.
aetsoil = 0.

!-------------------------
! Calculate actual canopy conductance, potential water scalar and actual evapotranspiration (aet) for each ipft

do pft = 1, npft

  if (present(pft)) then

    ! rootprop(1) = pftpar(1,pft) * 0.3
    ! rootprop(2) = pftpar(1,pft) * 0.45
    ! rootprop(3) = pftpar(1,pft) * 0.25
    ! rootprop(4) = (1.0 - pftpar(1,pft)) * 0.3
    ! rootprop(5) = (1.0 - pftpar(1,pft)) * 0.3
    ! rootprop(6) = (1.0 - pftpar(1,pft)) * 0.1

    rootprop(1) = rootfracl(pft,1)
    rootprop(2) = rootfracl(pft,2)
    rootprop(3) = rootfracl(pft,3)
    rootprop(4) = rootfracl(pft,4)
    rootprop(5) = rootfracl(pft,5)
    rootprop(6) = rootfracl(pft,6)

    ! Calculate effective supply function in root zone today, Eq. 24, Haxeltine & Prentice 1996

    wr = rootprop(1) * swf(1) + rootprop(2) * swf(2) + rootprop(3) * swf(3) + &
         rootprop(4) * swf(4) + rootprop(5) * swf(5) + rootprop(6) * swf(6)

    supply = emax * wr * ((rm_ind(pft) / sum(rm_ind)) ** (2./3.))

    ! Calculate actual evapotranspiration demand function and potential demand assuming full leaf cover
    ! Eqn 23, Haxeltine & Prentice 1996
    if (dphen(pft) > 0.) then

      demand = dpet * alpham * (1. - exp(-dgp(pft) * dphen(pft) / gm))

      !alternative formulation used by Gerten et al. 2004 (results in slightly lower demand for any given dgp)
      !demand = (1. - wet) * dpet(d) * alpham / (1. + gm / (dgp(d,pft) * dphen(d,pft)))

    else
      demand = 0.
    end if

    if (dphen(pft) == 1.) then
      demandpot = demand
    else
      demandpot = dpet * alpham * (1. - exp(-dgp(pft) / gm))
    end if

    !Calculate daily potential water scalar
    if (demandpot > 1.e-10) then                 !FLAG changed here from demand to demandpot to avoid flucutating leafout situation
      dwscal(pft) = min(supply / demandpot, 1.)
    else
      dwscal(pft) = 1.
    end if

    !--------------------------
    !if (ipft == 9) then
    !  write(stdout,'(a,4f9.4)')'waterbalance dwscal out',supply,demandpot,supply/demandpot,dwscal(ipft)
    !end if
    !--------------------------

    ! Calculate actual canopy conductance (gc), according to balance between supply and demand Eq. 25, Haxeltine & Prentice 1996

    if (supply >= demand) then
      dgc(pft) = dgp(pft) * dphen(pft)
    else if (dpet > 0.) then
      dgc(pft) = max(-gm * log(1. - supply / (dpet * alpham)),0.)
    else
      dgc(pft) = 0.
    end if

    !AET is smaller of supply and demand
    !Eqn 22, Haxeltine & Prentice 1996

    aet(pft) = min(supply,demand)

    daet = daet + aet(pft)

    !if (idx == 1) then
    !  write(*,*)ipft,aet(ipft),supply,demand
    !end if


    !Accumulate total AET
    !Eqns 28,29, Haxeltine & Prentice 1996

    if (wr == 0.) then
      beta(1) = 0.
      beta(2) = 0.
      beta(3) = 0.
      beta(4) = 0.
      beta(5) = 0.
      beta(6) = 0.
    else
      beta(1) = rootprop(1) * swf(1) / wr
      beta(2) = rootprop(2) * swf(2) / wr
      beta(3) = rootprop(3) * swf(3) / wr
      beta(4) = rootprop(4) * swf(4) / wr
      beta(5) = rootprop(5) * swf(5) / wr
      beta(6) = rootprop(6) * swf(6) / wr
    end if

    aetsoil(1) = aetsoil(1) + beta(1) * aet(pft)
    aetsoil(2) = aetsoil(2) + beta(2) * aet(pft)
    aetsoil(3) = aetsoil(3) + beta(3) * aet(pft)
    aetsoil(4) = aetsoil(4) + beta(4) * aet(pft)
    aetsoil(5) = aetsoil(5) + beta(5) * aet(pft)
    aetsoil(6) = aetsoil(6) + beta(6) * aet(pft)

    ! print *, pft, wr, supply, demand, swf

  end if

end do

!-------------------------
! Test
! Try invdiviodual PFT bucket? Do plant compete for water?

if (daet > dpet) then

  ! print *, ' daet > dpet'

  do pft = 1, npft

    aet(pft) = (aet(pft) / daet) * dpet

  end do

  daet = dpet

end if

!-------------------------

!-------------------------
!put this statement in to approximate some evaporation from the soil surface, even if there are no plants
!fpc_grid is used to estimate the bare ground fraction

! aetsoil(1) = aetsoil(1) + swf(1) * dpet * (1. - sum(fpc_grid))

end subroutine waterbalance

!---------------------------------------------------------------------

! subroutine soilwaterLPJ(awc,ksat,prec,aetsoil,swf,runoff,drainage) !,melt,mat20,idx)
!
! implicit none
!
! ! integer(i8) :: idx
!
! !arguments
!
! real(sp), intent(in),    dimension(2) :: awc       !soil available water content (mm or kg)
! real(sp), intent(in),    dimension(2) :: ksat      !saturated conductivity (mm h-1)
!
! ! real(sp), intent(in)                  :: melt      !daily total snowmelt (mm)
! real(sp), intent(in)                  :: prec      !daily total precipitation (mm)
!
! real(sp), intent(in),    dimension(2) :: aetsoil  !total water removed through evapotranspiration from each soil layer (mm)
!
! real(sp), intent(inout), dimension(2) :: swf       !water status in each soil layer (fraction of awc)
!
! real(sp), intent(out)                 :: runoff    !surface runoff (mm)
! real(sp), intent(out)                 :: drainage  !groundwater drainage (mm)
!
! ! real(sp), intent(in) :: mat20  !20 year running mean annual temperature
!
! !parameters
!
! integer,  parameter :: ts = 24
! real(sp), parameter :: dt = 1. / real(ts)
! real(sp), parameter :: a  = 1.e-5  !minimum unsaturated conductivity (mm hr-1)
!
! !local variables
!
! integer  :: t         !counters
! real(sp) :: qsurf     !surface flux of water (mm h-1)
! real(sp) :: perc      !percolation of soil water (mm) between layers
! real(sp) :: reversef  !reverse flow: if lower layer is saturated, water backs up into the upper layer
!
! real(sp) :: infil
!
!
! real(sp), dimension(2) :: whcmm   !water holding capacity of the soil in each layer (mm)
! real(sp), dimension(2) :: kmmh    !saturated conductivity of the soil in each layer (mm hr-1)
! real(sp), dimension(2) :: swater  !soil water content (mm)
! real(sp), dimension(2) :: aet     !evapotranspiration (mm h-1)
!
! !------------------------------
!
! infil    = 0.
! perc     = 0.
! drainage = 0.
! runoff   = 0.
! reversef = 0.
!
! whcmm = awc   !(/  84.56579590,  300.7981567   /)
! kmmh  = ksat  !(/   7.46948719,    5.689514160 /)
!
! !instantaneous hydraulic conductivity is based on the power4 scaling to wetness fraction (soil at FC behaves as saturated flow)
!
! ! qsurf = dt * (melt + prec)   !mm hr-1
! qsurf = dt * prec   !mm hr-1
! aet   = dt * aetsoil
!
! do t = 1, ts  !hourly loop
!
!   kmmh = ksat * swf**4 - a*swf**4 + a  !unsaturated conductivity (mm hr-1)
!
!   swater = swf * whcmm  !total water mass (mm)
!
!   !drainage out the bottom
!
!   ! if (mat20 <= 0.) then    !very simple parameterization of permafrost impedance of drainage
!   !   drainage = 0.
!   ! else
!   !   drainage = min(swater(2),kmmh(2))
!   ! end if
!
!   drainage = min(swater(2),kmmh(2))
!
!   !percolation between layers
!
!   perc = min(swater(1),kmmh(1))
!
!   !infiltration is a function of the available pore volume and the conductivity
!
!   infil = qsurf
!
!   swater(1) = max(swater(1) + infil - aet(1) - perc, 0.)
!   swater(2) = max(swater(2) + perc  - aet(2) - drainage, 0.)
!
!   if (swater(2) > whcmm(2)) then
!     reversef  = swater(2) - whcmm(2)
!     swater(1) = swater(1) + reversef
!     swater(2) = whcmm(2)
!   end if
!
!   if (swater(1) > whcmm(1)) then
!     runoff = runoff + swater(1) - whcmm(1)
!     swater(1) = whcmm(1)
!   end if
!
!   !update swf
!
!   where (whcmm > 0.)
!     swf = swater / whcmm
!   elsewhere
!     swf = 0.
!   end where
!
! end do
!
! !if (idx == 1) then
! !  write(*,*)melt,prec,runoff,drainage,aetsoil
! !end if
!
!
! if (swf(2) <= 0.) then
!
!   !write(stdout,'(5f10.4)')prec,melt,qsurf,infil,runoff
!   !write(stdout,'(5f10.4)')whcmm(1),kmmh(1),swater(1),swf(1),perc
!   !write(stdout,'(6f10.4)')whcmm(2),kmmh(2),swater(2),swf(2),drainage,reversef
!
! !  stop
! end if
!
! end subroutine soilwaterLPJ

!---------------------------------------------------------------------

end module waterbalancemod
