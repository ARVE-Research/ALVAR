module waterbalancemod

! Water balance module copied from LPJ-LMFire (Feb 2021 by Leo O Lai)
! Include alterantive simplesoilLPJ() subroutine for two layer soil model

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec
use statevarsmod,  only : cnt

implicit none

real(sp), allocatable, dimension(:,:) :: soilpar
real(sp), allocatable, dimension(:,:) :: w

integer(i4), parameter :: nl_LPJ = 2

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine simplesoilLPJ(grid)

use statevarsmod, only : soilvars
use soilstatemod, only : calctheta,fbulk,fKsat

implicit none

integer(i4), intent(in) :: grid

! Parameters
real(sp), parameter :: ombd = 0.224  !bulk density of organic matter (g cm-3)
real(sp), parameter :: omcf = 1.724  !conversion factor from organic carbon to organic matter

integer :: l
integer :: it

integer :: soiltype

real(sp) :: sand
real(sp) :: clay
real(sp) :: OM    !(mass %)
real(sp) :: cfvo

real(sp) :: silt
real(sp) :: bulk
real(sp) :: Tsat
real(sp) :: T33
real(sp) :: T1500

real(sp), allocatable, dimension(:) :: dz
real(sp), allocatable, dimension(:) :: zpos
real(sp), allocatable, dimension(:) :: OrgM  !(g m-2)
real(sp), allocatable, dimension(:) :: whc
real(sp), allocatable, dimension(:) :: Ksat

real(sp) :: Csoil     !(g m-2)
real(sp) :: orgC      !(mass %)
real(sp) :: soilmass  !(g m-3)
real(sp) :: blk0      !(g m-3)
real(sp) :: dOM       !change in organic matter (g m-2)
real(sp) :: dzOM      !interlayer transport of SOM (g m-2)
real(sp) :: dzx       !excess change in top layer thickness

logical, allocatable, dimension(:) :: valid

!-------------------------

if (grid == 1) allocate(soilpar(cnt(1),8))
if (grid == 1) allocate(w(cnt(1),2))

allocate(dz(nl_LPJ))
allocate(zpos(nl_LPJ))
allocate(OrgM(nl_LPJ))
allocate(whc(nl_LPJ))
allocate(Ksat(nl_LPJ))
allocate(valid(nl_LPJ))

valid = .true.

zpos = [0.15,0.45]
dz   = 0.3

OrgM = 0.

!-------------------------

do l = 1, nl_LPJ

  if (l == 1) then
    sand  = sum(soilvars(grid)%sand(1:3)) / 3.
    clay  = sum(soilvars(grid)%clay(1:3)) / 3.
    cfvo  = sum(soilvars(grid)%cfvo(1:3)) / 3.
    OM    = sum(soilvars(grid)%OrgM(1:3)) / 3.
  else if (l == 2) then
    sand  = soilvars(grid)%sand(4)
    clay  = soilvars(grid)%clay(4)
    cfvo  = soilvars(grid)%cfvo(4)
    OM    = soilvars(grid)%OrgM(4)
  end if

  silt = 100. - (sand + clay)

!   write(0,*)'SOIL',l,sand,silt,clay,OM

  if (sand < 0.) then
    valid(l) = .false.
    cycle
  end if

  Csoil = OrgM(l) / omcf  !layer

  orgC = Csoil

  !write(stdout,*)'lyr ',l,' tile',i
  !write(stdout,*)'dz  ',dz(l)
  !write(stdout,*)'Csol',Csoil

  !because bulk density depends strongly on organic matter content and
  !weakly on wilting point water content, we guess an initial value and
  !iterate to a stable solution

  T1500 = 0.1

  bulk = fbulk(orgC,T1500*100.,clay,zpos(l),silt)

  it = 1

  do

    !convert SOM to mass fraction and calculate the difference in OM content

    soilmass = bulk * 1.e6 * dz(l)     !(g) (m-2)  note: g cm-3 * (m*100=cm) 100cm * 100cm = g

    orgC = max(100. * Csoil / soilmass,0.)  !(mass %)

    !recalculate bulk

    blk0 = fbulk(orgC,T1500*100.,clay,zpos(l),silt)  !units (g cm-3)

    !calculate wilting point, field capacity, and saturation, needs input in fractions not percent

    OM = orgC * omcf             !(mass %)

    if (OM >= 30.) soiltype = 3  !humic soil

    call calctheta(sand/100.,clay/100.,OM/100.,blk0,Tsat,T33,T1500,soiltype)

    !recalculate bulk

    bulk = fbulk(orgC,T1500*100.,clay,zpos(l),silt)  !units (g cm-3)

    if (abs(bulk - blk0) < 0.001 .or. it > 50) exit

    blk0 = bulk

    it = it + 1

  end do

  !with the final value for bulk density, recalculate porosity

  call calctheta(sand/100.,clay/100.,OM/100.,bulk,Tsat,T33,T1500)

  !update layer-integrated WHC

  whc(l) = 1000. * dz(l) * (T33 - T1500)  !mm

  !calculate saturated conductivity

  Ksat(l) = fKsat(Tsat,T33,T1500)

  !write(stdout,*)'layer',l
  !write(stdout,*)'sand',sand
  !write(stdout,*)'silt',silt
  !write(stdout,*)'clay',clay
  !write(stdout,*)'OM  ',OM
  !write(stdout,*)'bulk',bulk

  !write(stdout,*)'Tsat',Tsat
  !write(stdout,*)'T33 ',T33
  !write(stdout,*)'Twp ',T1500
  !write(stdout,*)'Ksat',Ksat(l)

  soilpar(grid,8) = bulk

end do  !layers

!-------------------------

where (.not. valid)  !assign a typical rock value for non-soil layers
  whc  = 0.03
  Ksat = 1.e-3
  valid = .true.
end where

!saturated conductivity (mm h-1)
soilpar(grid,1) = Ksat(1)
soilpar(grid,2) = sum(Ksat(2:nl_LPJ))/count(valid(2:nl_LPJ))

!water holding capacity (mm)
soilpar(grid,3) = whc(1)
soilpar(grid,4) = sum(whc(2:nl_LPJ),mask=valid(2:nl_LPJ)) + whc(nl_LPJ) / dz(nl_LPJ) * 2. !mm/m layer add 2 more meters like LPJ2

soilpar(grid,5) = 0.2    !thermal diffusivity (mm2/s) at wilting point (0% WHC)
soilpar(grid,6) = 0.650  !thermal diffusivity (mm2/s) at 15% WHC
soilpar(grid,7) = 0.4    !thermal diffusivity at field capacity (100% WHC)

end subroutine simplesoilLPJ

!---------------------------------------------------------------------

subroutine waterbalanceLPJ(grid,d,present,dgp,dpet,dphen,  &
                        dgc,dprec,drunoff_drain,  &
                        drunoff_surf,dwscal,daet,fpc_grid) !,w,dmelt,ksat,awc,rootprop,mat20)

use pftparmod, only : npft,pftpar

implicit none

! Arguments
integer,  intent(in) :: grid
integer,  intent(in) :: d
! real(sp), intent(in) :: mat20  !20 year running mean annual temperature
real(sp), intent(in) :: dpet
! real(sp), intent(in) :: dmelt
real(sp), intent(in) :: dprec
! real(sp), intent(in) :: ksat  !(mm hr-1)
! real(sp), intent(in) :: awc
real(sp), dimension(:), intent(in) :: fpc_grid
logical, dimension(:),  intent(in) :: present
real(sp), intent(in) :: dphen  !phenological state (0=no leaves, 1=full canopy)
real(sp), dimension(:), intent(in) :: dgp    !non water stressed canopy conductance (?)

! real(sp), intent(inout), dimension(:) :: w

real(sp), intent(out) :: daet
real(sp), intent(out) :: drunoff_surf
real(sp), intent(out) :: drunoff_drain
real(sp), dimension(:), intent(out) :: dwscal
real(sp), dimension(:), intent(out) :: dgc    !actual canopy conductance (?)

! Parameters
real(sp), parameter :: emax   = 5.   !maximum daily transpiration rate (mm/day)
real(sp), parameter :: alpham = 1.4  !maximum Priestly-Taylor coefficient
real(sp), parameter :: gm     = 5.   !scaling conductance (mm s-1)?

! Local variables
real(sp), dimension(2) :: rootprop
real(sp), dimension(2) :: ksat
real(sp), dimension(2) :: awc
integer :: pft

real(sp) :: supply
real(sp) :: demand
real(sp) :: demandpot
real(sp) :: wr

real(sp), dimension(2) :: aettotal
real(sp), dimension(2) :: beta

real(sp), dimension(npft) :: aet

!-------------------------

daet = 0.
aettotal = 0.
drunoff_surf = 0.
drunoff_drain = 0.

ksat(1) = soilpar(grid,1)
ksat(2) = soilpar(grid,2)
awc(1)  = soilpar(grid,3)
awc(2)  = soilpar(grid,4)

!-------------------------
!Calculate actual canopy conductance, potential water scalar and actual evapotranspiration (aet) for each ipft

do pft = 1, npft

  if (present(pft)) then

    rootprop(1) = pftpar(1,pft)
    rootprop(2) = 1.0 - rootprop(1)

    !Calculate effective supply function in root zone today, Eq. 24, Haxeltine & Prentice 1996

    wr = rootprop(1) * w(grid,1) + rootprop(2) * w(grid,2)

    supply = emax * wr

    !Calculate actual evapotranspiration demand function and potential demand assuming full leaf cover
    !Eqn 23, Haxeltine & Prentice 1996

    if (dphen > 0.) then

      demand = dpet * alpham * (1. - exp(-dgp(pft) * dphen / gm))

      !alternative formulation used by Gerten et al. 2004 (results in slightly lower demand for any given dgp)
      !demand = (1. - wet) * dpet(d) * alpham / (1. + gm / (dgp(d,pft) * dphen(d,pft)))

    else
      demand = 0.
    end if

    if (dphen == 1.) then
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
      dgc(pft) = dgp(pft) * dphen
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
    else
      beta(1) = rootprop(1) * w(grid,1) / wr
      beta(2) = rootprop(2) * w(grid,2) / wr
    end if

    aettotal(1) = aettotal(1) + beta(1) * aet(pft)
    aettotal(2) = aettotal(2) + beta(2) * aet(pft)

  end if

end do

!-------------------------
!put this statement in to approximate some evaporation from the soil surface, even if there are no plants
!fpc_grid is used to estimate the bare ground fraction

aettotal(1) = aettotal(1) + w(grid,1) * dpet * (1. - sum(fpc_grid))

call soilwaterLPJ(awc,ksat,dprec,aettotal,w(grid,:),drunoff_surf,drunoff_drain)

!write(90,*)year,d,dpet(d),sum(aettotal),dprec(d)

end subroutine waterbalanceLPJ

!---------------------------------------------------------------------

subroutine soilwaterLPJ(awc,ksat,prec,aettotal,swf,runoff,drainage) !,melt,mat20,idx)

implicit none

! integer(i8) :: idx

!arguments

real(sp), intent(in),    dimension(2) :: awc       !soil available water content (mm or kg)
real(sp), intent(in),    dimension(2) :: ksat      !saturated conductivity (mm h-1)

! real(sp), intent(in)                  :: melt      !daily total snowmelt (mm)
real(sp), intent(in)                  :: prec      !daily total precipitation (mm)

real(sp), intent(in),    dimension(2) :: aettotal  !total water removed through evapotranspiration from each soil layer (mm)

real(sp), intent(inout), dimension(2) :: swf       !water status in each soil layer (fraction of awc)

real(sp), intent(out)                 :: runoff    !surface runoff (mm)
real(sp), intent(out)                 :: drainage  !groundwater drainage (mm)

! real(sp), intent(in) :: mat20  !20 year running mean annual temperature

!parameters

integer,  parameter :: ts = 24
real(sp), parameter :: dt = 1. / real(ts)
real(sp), parameter :: a  = 1.e-5  !minimum unsaturated conductivity (mm hr-1)

!local variables

integer  :: t         !counters
real(sp) :: qsurf     !surface flux of water (mm h-1)
real(sp) :: perc      !percolation of soil water (mm) between layers
real(sp) :: reversef  !reverse flow: if lower layer is saturated, water backs up into the upper layer

real(sp) :: infil


real(sp), dimension(2) :: whcmm   !water holding capacity of the soil in each layer (mm)
real(sp), dimension(2) :: kmmh    !saturated conductivity of the soil in each layer (mm hr-1)
real(sp), dimension(2) :: swater  !soil water content (mm)
real(sp), dimension(2) :: aet     !evapotranspiration (mm h-1)

!------------------------------

infil    = 0.
perc     = 0.
drainage = 0.
runoff   = 0.
reversef = 0.

whcmm = awc   !(/  84.56579590,  300.7981567   /)
kmmh  = ksat  !(/   7.46948719,    5.689514160 /)

!instantaneous hydraulic conductivity is based on the power4 scaling to wetness fraction (soil at FC behaves as saturated flow)

! qsurf = dt * (melt + prec)   !mm hr-1
qsurf = dt * prec   !mm hr-1
aet   = dt * aettotal

do t = 1, ts  !hourly loop

  kmmh = ksat * swf**4 - a*swf**4 + a  !unsaturated conductivity (mm hr-1)

  swater = swf * whcmm  !total water mass (mm)

  !drainage out the bottom

  ! if (mat20 <= 0.) then    !very simple parameterization of permafrost impedance of drainage
  !   drainage = 0.
  ! else
  !   drainage = min(swater(2),kmmh(2))
  ! end if

  drainage = min(swater(2),kmmh(2))

  !percolation between layers

  perc = min(swater(1),kmmh(1))

  !infiltration is a function of the available pore volume and the conductivity

  infil = qsurf

  swater(1) = max(swater(1) + infil - aet(1) - perc, 0.)
  swater(2) = max(swater(2) + perc  - aet(2) - drainage, 0.)

  if (swater(2) > whcmm(2)) then
    reversef  = swater(2) - whcmm(2)
    swater(1) = swater(1) + reversef
    swater(2) = whcmm(2)
  end if

  if (swater(1) > whcmm(1)) then
    runoff = runoff + swater(1) - whcmm(1)
    swater(1) = whcmm(1)
  end if

  !update swf

  where (whcmm > 0.)
    swf = swater / whcmm
  elsewhere
    swf = 0.
  end where

end do

!if (idx == 1) then
!  write(*,*)melt,prec,runoff,drainage,aettotal
!end if


if (swf(2) <= 0.) then

  !write(stdout,'(5f10.4)')prec,melt,qsurf,infil,runoff
  !write(stdout,'(5f10.4)')whcmm(1),kmmh(1),swater(1),swf(1),perc
  !write(stdout,'(6f10.4)')whcmm(2),kmmh(2),swater(2),swf(2),drainage,reversef

!  stop
end if

end subroutine soilwaterLPJ

!---------------------------------------------------------------------

end module waterbalancemod
