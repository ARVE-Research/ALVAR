module simplesoilmod

use parametersmod, only : i2,i4,sp,dp,missing_sp

implicit none

public :: simplesoil

private :: fbulk
private :: calctheta
private :: fKsat
private :: fTsat

!parameters
integer(i4), parameter :: nl = 6
real(sp), parameter :: ombd = 0.224  !bulk density of organic matter (g cm-3)
real(sp), parameter :: omcf = 1.724  !conversion factor from organic carbon to organic matter

! Soil layer depth from soil input data file
real(sp), dimension(nl) :: zpos

!---------------------------------------------------------------------

type soildata

  real(sp), dimension(nl) :: sand
  real(sp), dimension(nl) :: clay
  real(sp), dimension(nl) :: cfvo
  real(sp), dimension(nl) :: OrgM
  real(sp), dimension(nl) :: bulk
  real(sp), dimension(nl) :: whc
  real(sp), dimension(nl) :: Ksat

end type soildata

type(soildata), target, allocatable, dimension(:) :: soilvars

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine simplesoil(grid)

! use mpistatevarsmod, only : soildata
! use pedotransfermod, only : fbulk,calctheta,fKsat,ombd,omcf

implicit none

integer(i4), intent(in) :: grid

! type(soildata), intent(inout) :: soil  !state variables sent back out with MPI
! real(sp), dimension(:), intent(out) :: soilpar

integer(i4) :: l
integer(i4) :: it

integer(i4) :: soiltype

real(sp) :: sand
real(sp) :: clay
real(sp) :: OM    !(mass %)

real(sp) :: silt
real(sp) :: bulk
real(sp) :: Tsat
real(sp) :: T33
real(sp) :: T1500

real(sp), allocatable, dimension(:) :: dz
logical,  allocatable, dimension(:) :: valid

! real(sp), pointer, dimension(:) :: zpos
real(sp), pointer, dimension(:) :: OrgM  !(g m-2)
real(sp), pointer, dimension(:) :: whc
real(sp), pointer, dimension(:) :: Ksat

real(sp) :: Csoil     !(g m-2)
real(sp) :: orgC      !(mass %)
real(sp) :: soilmass  !(g m-3)
real(sp) :: blk0      !(g m-3)
real(sp) :: dOM       !change in organic matter (g m-2)
real(sp) :: dzOM      !interlayer transport of SOM (g m-2)
real(sp) :: dzx       !excess change in top layer thickness

!
! !----------
!
! nl = size(soil%sand)
!
! ! write(0,*)'SOIL NLAYERS',nl

OrgM => soilvars(grid)%OrgM
whc  => soilvars(grid)%whc
Ksat => soilvars(grid)%Ksat

allocate(dz(nl))
allocate(valid(nl))
! allocate(zpos(nl))
! allocate(OrgM(nl))
! allocate(whc(nl))
! allocate(Ksat(nl))


valid = .true.

! Midpoint of soil layer zpos = [-2.5, -10, -22.5, -45, -80, -150] cm
! dz = soil layer depth in meter
dz = [0.05, 0.1, 0.15, 0.3, 0.4, 1.0]

OrgM = 1.

soiltype = 1

layerloop : do l = 1, nl

  sand  = soilvars(grid)%sand(l)
  clay  = soilvars(grid)%clay(l)
  OM    = soilvars(grid)%OrgM(l)

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

  ! if (it > 1) print *, it

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

  soilvars(grid)%bulk(l) = bulk

end do layerloop


where (.not. valid)  !assign a typical rock value for non-soil layers
  ! whc  = 0.03
  ! Ksat = 1.e-3
  ! valid = .true.

  soilvars(grid)%whc  = missing_sp
  soilvars(grid)%Ksat = missing_sp
  soilvars(grid)%bulk = missing_sp
  valid = .true.
end where


! !saturated conductivity (mm h-1)
! soilpar(1) = Ksat(1)
! soilpar(2) = sum(Ksat(2:nl))/count(valid(2:nl))
!
! !water holding capacity (mm)
! soilpar(3) = whc(1)
! soilpar(4) = sum(whc(2:nl),mask=valid(2:nl)) + whc(nl) / dz(nl) * 2. !mm/m layer add 2 more meters like LPJ2
!
! soilpar(5) = 0.2    !thermal diffusivity (mm2/s) at wilting point (0% WHC)
! soilpar(6) = 0.650  !thermal diffusivity (mm2/s) at 15% WHC
! soilpar(7) = 0.4    !thermal diffusivity at field capacity (100% WHC)


end subroutine simplesoil

!---------------------------------------------------------------------

subroutine calctheta(sand,clay,OM,bulk,Tsat,T33,T1500,soiltype)

! equations from JAP Pollacco, 2008. Can. J. Soil Sci. 88: 761-774
! subroutine copied from LPJ-LMFire (Leo Lai June 2021)

implicit none

!arguments

real(sp), intent(in)  :: sand               !mass fraction of sand (g g-1) (0-1)
real(sp), intent(in)  :: clay               !mass fraction of clay (g g-1) (0-1)
real(sp), intent(in)  :: OM                 !mass fraction of organic matter (g g-1) (0-1) (multiply by 1.724 if input data is in terms of Carbon)
real(sp), intent(in)  :: bulk               !total soil bulk density (g cm-3)
integer,  intent(in), optional :: soiltype  !soil type code 1=typical, 2=tropical, 3=humic, 4=vitric

real(sp), intent(out) :: Tsat               !volumetric saturated water content (m3 m-3)
real(sp), intent(out) :: T33                !volumetric water content at field capacity (Psi = -33 kPa)   (m3 m-3)
real(sp), intent(out) :: T1500              !volumetric water content at wilting point  (Psi = -1500 kPa) (m3 m-3)

!parameters

type fitpars
  real(sp), dimension(2) :: Pmax
  real(sp), dimension(2) :: Pmin
  real(sp), dimension(2) :: Psand
  real(sp), dimension(2) :: Pclay
  real(sp), dimension(2) :: Pp
end type fitpars

type(fitpars) :: typical
type(fitpars) :: tropical
type(fitpars) :: humic
type(fitpars) :: vitric

type(fitpars) :: pars

!local variables

real(sp) :: sand_om   !sand fraction corrected for organic matter (0-1) (sand_om = (1-OM)*sand)
real(sp) :: clay_om

real(sp), dimension(2) :: Pmin
real(sp), dimension(2) :: Pmax
real(sp), dimension(2) :: Pclay
real(sp), dimension(2) :: Psand
real(sp), dimension(2) :: Pp

real(sp) :: Wsat    !saturated gravimetric soil moisture content (g g-1)

real(sp), dimension(2) :: W  !gravimetric soil moisture content (g g-1)

!---------------------------------------------------------------------------
!optimal fitting parameters for W33 and W1500, from Table 7 in Pollacco 2008

!                    FC     WP
typical%Pmax   = [ 0.953, 0.895 ]
typical%Pmin   = [ 0.608, 0.165 ]
typical%Psand  = [ 0.215, 0.000 ]
typical%Pclay  = [ 0.914, 0.759 ]
typical%Pp     = [-0.102, 1.468 ]

!                    FC     WP
tropical%Pmax  = [ 1.000, 0.891 ]
tropical%Pmin  = [ 0.781, 0.197 ]
tropical%Psand = [ 0.338, 0.000 ]
tropical%Pclay = [ 2.104, 0.521 ]
tropical%Pp    = [-2.009, 0.767 ]

!                    FC     WP
humic%Pmax     = [ 0.685, 0.551 ]
humic%Pmin     = [ 0.654, 0.190 ]
humic%Psand    = [ 0.217, 0.000 ]
humic%Pclay    = [ 3.010, 0.372 ]
humic%Pp       = [-1.810, 0.273 ]

!                    FC     WP
vitric%Pmax    = [ 1.000, 1.000 ]
vitric%Pmin    = [ 0.371, 0.094 ]
vitric%Psand   = [ 0.187, 0.000 ]
vitric%Pclay   = [ 0.563, 0.757 ]
vitric%Pp      = [-0.030, 0.616 ]

!---------------------------------------------------------------------------

if (present(soiltype)) then
  select case(soiltype)
  case default
    pars = typical
  case(2)
    pars = tropical
  case(3)
    pars = humic
  case(4)
    pars = vitric
  end select
else
  pars = typical
end if

Pmin  = pars%Pmin
Pmax  = pars%Pmax
Psand = pars%Psand
Pclay = pars%Pclay
Pp    = pars%Pp

!-----

sand_om = (1. - OM) * sand
clay_om = (1. - OM) * clay

Tsat = fTsat(bulk,OM,sand_om)

if (bulk > 0.) then

  Wsat = Tsat / bulk

  !Pollacco PTF Model 4, eqn. 7a
  W = Wsat * (Pmin + (Pmax - Pmin) * clay_om**(Pclay + Pp * Wsat**2)) * exp(-(Psand * sand_om**3) / Wsat)

else

  Wsat = 0.
  W    = 0.

end if

T33   = W(1) * bulk   !eqn. 1, NB pw (water density assumed = 1)
T1500 = W(2) * bulk


end subroutine calctheta

!---------------------------------------------------------------------------

function fTsat(bulk,OM,sand_om)

!equations from JAP Pollacco, 2008. Can. J. Soil Sci. 88: 761-774

implicit none

real(sp) :: fTsat

real(sp), intent(in) :: bulk     !bulk density (g cm-3)
real(sp), intent(in) :: OM       !organic matter mass fraction (0-1)
real(sp), intent(in) :: sand_om  !sand fraction corrected for organic matter (0-1) (sand_om = (1-OM)*sand)

real(sp) :: pmin   !mineral density (g cm-3)
real(sp) :: pom    !particle density of organic matter (g cm-3)
real(sp) :: pp     !soil particle density (g cm-3)

!-----

pmin  = 2.69 - 0.03 * sand_om  !pg 763

pom   = 1.127 + 0.373 * OM     !pg. 763

pp    = 1. / ((OM / pom) + ((1. - OM) / pmin))  !eqn. 3

fTsat = 1. - bulk / pp   !eqn. 2a

end function fTsat

!---------------------------------------------------------------------------

function fbulk(orgC,W15,clay,depth,silt)

!equations from Heuscher et al., (2005) Soil Sci. Soc. Am. J. 69

real(sp) :: fbulk  !bulk density (g cm-3)

real(sp), intent(in) :: orgC   !mass %
real(sp), intent(in) :: W15    !water content at -15 bar (-1500 kPa) (mass %)
real(sp), intent(in) :: clay   !mass %
real(sp), intent(in) :: silt   !mass %
real(sp), intent(in) :: depth  !horizon depth (cm)

real(sp), parameter :: incp  =  1.685    !intercept coefficient
real(sp), parameter :: pOC   = -0.198    !organic carbon coefficient
real(sp), parameter :: pwc   = -0.0133   !wilting point water content coefficient
real(sp), parameter :: pclay =  0.0079   !clay coef
real(sp), parameter :: pdep  =  0.00014  !depth coef
real(sp), parameter :: psilt = -0.0007   !silt coef

fbulk = incp + pOC * sqrt(orgC) + pwc * W15 + pclay * clay + pdep * depth + psilt * silt

end function fbulk

!---------------------------------------------------------------------------

function fKsat(Tsat,T33,T1500)

!equations from Saxton & Rawls (2006) Soil Sci. Soc. Am. J. 70

real(sp) :: fKsat  !(mm h-1)

real(sp), intent(in) :: Tsat
real(sp), intent(in) :: T33
real(sp), intent(in) :: T1500

real(sp), parameter :: l1500 = log(1500.)
real(sp), parameter :: l33   = log(33.)
real(sp), parameter :: num   = l1500 - l33

real(sp) :: B
real(sp) :: lambda

B = num / (log(T33) - log(T1500))

lambda = 1. / B

fKsat = 1930. * (Tsat - T33)**(3.-lambda)

end function fKsat

!---------------------------------------------------------------------------


end module simplesoilmod
