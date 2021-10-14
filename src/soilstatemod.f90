module soilstatemod

! Module copied from LPJ-LMFire and ARVE-DGVM

use parametersmod, only : i2,i4,sp,dp,missing_sp,Tfreeze

implicit none

public :: soilprep

private :: fbulk
private :: fRock
private :: calctheta
private :: fKsat
private :: fTsat
private :: fPsat
private :: fBexp

! Number of soil layers
integer(i4), parameter               :: nl = 6

! Soil layer thickness
real(sp),    parameter, dimension(6) :: dz = [0.05, 0.1, 0.15, 0.3, 0.4, 1.0]
real(sp),    parameter, dimension(6) :: dzmm = dz * 1.e3

! Midpoint z position (depth) of soil layer
real(sp),    parameter, dimension(6) :: zpos = [0.025, 0.1, 0.225, 0.45, 0.8, 1.5]
real(sp),    parameter, dimension(6) :: zposmm = zpos * 1.e3

! Soil layer interface z position (depth), positive downwards (dim = nl+1, including surface and column bottom interface)
real(sp),    parameter, dimension(7) :: zipos = [0., 0.05, 0.15, 0.3, 0.6, 1.0, 2.0]
real(sp),    parameter, dimension(7) :: ziposmm = zipos * 1.e3


!-------------------------

! Soil properties parameters (copied from LPJ-LMFire, Leo Lai Jun 2021)
real(sp), parameter :: ombd = 0.224  !bulk density of organic matter (g cm-3)
real(sp), parameter :: omcf = 1.724  !conversion factor from organic carbon to organic matter

! Soil properties parameters (copied from ARVE-DGVM, Leo Lai Jul 2021)
real(dp), parameter :: OMorgC = 0.5800464               ! g org C =sorg/ 1.724 (conversion factor from Nelson & Sommers 1996) !FLAG THIS IS FOR SOILS!
real(dp), parameter :: peatlim= 25._dp                  ! percent organic matter in soil where the soil is treated as peat
real(dp), parameter :: soilbulk = 2650._dp              ! soil bulk density (kg m-3)
real(dp), parameter :: omd = 1.3                        ! (g cm-3) particle density of org matter (from Hillel 1982)
real(dp), parameter :: sandd = 2.128_dp                 ! (g cm-3) particle density of sand(0.8 x quartz (2.66))
real(dp), parameter :: siltd = 2.35_dp                  ! (average soil density...   FLAG!!!must look up a better value for silt!!!
real(dp), parameter :: clayd = 2.65_dp                  ! (g cm-3) particle density of other minerals
real(dp), parameter :: rockd = 2.70_dp                  ! (g cm-3) particle density of coarse frags
real(dp), parameter :: sandquartz = 0.65                ! fraction of quartz in sand on average (geology.uprm.edu/Morelock/terrigenous.html)
real(dp), parameter :: rockquartz = 0.25                ! fraction of quartz in rock on average (estimate)
real(dp), parameter :: bedrock_tc = 3.7                 ! bedrock thermal conductivity (W m-1 K-1)
real(dp), parameter :: Corg =  2.51e6                   ! (1e6)Heat capacity of soil organic matter (typical) (J m-3 K-1)
real(dp), parameter :: Csand = 2.128e6                  ! (1e6) Volumetric heat capacity of sand (J m-3 K-1)
real(dp), parameter :: Csilt = 2.439e6                  ! (1e6) Volumetric heat capacity of silt (J m-3 K-1) from Ren et al SSSAJ 67 (6) 2003
real(dp), parameter :: Cclay = 2.8385e6                 ! (1e6) Volumetric heat capacity of clay (J m-3 K-1)
real(dp), parameter :: Crock = 2.01e6                   ! (1e6) Volumetric heat capacity of quartz,other minerals (J m-3 K-1)
real(dp), parameter :: Kom = 0.25_dp                    ! thermal cond of organic matter (W m-1 K-1)
real(dp), parameter :: Kom_dry = 0.05d0                 ! thermal cond of dry organic matter (W m-1 K-1) (Farouki, 1981)
real(dp), parameter :: Kquartz = 8.0_dp                 ! thermal cond. of quartz (W m-1 K-1) (applicable across more temps than previous value)
real(dp), parameter :: Kmineral = 2.5_dp                ! thermal cond. of other minerals (clay) (W m-1 K-1)
real(dp), parameter :: pliq  = 1000._dp                 ! density of water (kg m-3)

real(dp), parameter :: a = 0.053_dp                     ! From Balland and Arp for dry soil thermal conductivity calc
real(dp), parameter :: alpha = 0.24_dp                  ! +/- 0.04 constant for Knum calc.from Balland and Arp
real(dp), parameter :: beta = 18.3_dp                   ! +/- 1.1 constant for Knum calc.from Balland and Arp

real(dp), parameter :: Kliq = 0.57_dp                   ! Thermal conductivity of liquid water (typical) (W m-1 K-1)
real(dp), parameter :: Kice = 2.2_dp                    ! Thermal conductivity of ice (typical) (W m-1 K-1)
real(dp), parameter :: Kair = 0.0243_dp                 ! Thermal conductivity of dry air (W m-1 K-1)

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine soilprep(grid)

! Subroutine to calculate soil physical properties from sand, clay and cfvo soildata input (soil grid)
! Code copied (and combined) from LPJ-LMFire (simplesoilmod.f90) and ARVE-DGVM (init_model_state.f90)
! by Leo O Lai (Jul 2021)

use metvarsmod, only : soilvars

implicit none

integer(i4), intent(in) :: grid

! type(soildata), intent(inout) :: soil  !state variables sent back out with MPI
! real(sp), dimension(:), intent(out) :: soilpar

integer(i4) :: l
integer(i4) :: it

integer(i4) :: soiltype

real(sp) :: sand
real(sp) :: clay
real(sp) :: cfvo
real(sp) :: OM    !(mass %)

real(sp) :: silt
! real(sp) :: bulk1
! real(sp) :: Tsat1
real(sp) :: T33       ! Water content at field capacity
real(sp) :: T1500     ! Wilting point soil content

logical, allocatable, dimension(:) :: valid

! real(sp), pointer, dimension(:) :: zpos
logical, pointer :: validcell

real(sp), pointer, dimension(:) :: bulk
real(sp), pointer, dimension(:) :: rock  ! Course fragment content by mass (percent)
real(sp), pointer, dimension(:) :: OrgM  !(g m-2)

real(sp), pointer, dimension(:) :: Vsand
real(sp), pointer, dimension(:) :: Vclay
real(sp), pointer, dimension(:) :: Vsilt
real(sp), pointer, dimension(:) :: VOrgM

real(sp), pointer, dimension(:) :: whc
real(sp), pointer, dimension(:) :: Ksat
real(sp), pointer, dimension(:) :: Tsat
real(sp), pointer, dimension(:) :: Tfield
real(sp), pointer, dimension(:) :: Twilt
real(sp), pointer, dimension(:) :: Psat
real(sp), pointer, dimension(:) :: Bexp

real(sp), pointer, dimension(:) :: Csolid
real(sp), pointer, dimension(:) :: Ksolid
real(sp), pointer, dimension(:) :: Kdry

real(sp), pointer, dimension(:) :: Wliq
real(sp), pointer, dimension(:) :: Wice
real(sp), pointer, dimension(:) :: Tsoil
real(sp), pointer, dimension(:) :: Tsoiln
real(sp), pointer, dimension(:) :: Tliq
real(sp), pointer, dimension(:) :: Tice
real(sp), pointer, dimension(:) :: Psi

!local variables
real(sp) :: densp
real(sp) :: Vtotal
real(sp), dimension(1:nl) :: Kdrysolid    !soil dry mineral thermal conductivity at layer midpoint (W m-1 K-1)
real(sp) :: interm                        !temporary variable
real(sp) :: Vquartz                       !volume fraction quartz

real(sp) :: Csoil     !(g m-2)
real(sp) :: orgC      !(mass %)
real(sp) :: soilmass  !(g m-3)
real(sp) :: blk0      !(g m-3)
real(sp) :: dOM       !change in organic matter (g m-2)
real(sp) :: dzOM      !interlayer transport of SOM (g m-2)
real(sp) :: dzx       !excess change in top layer thickness

!-------------------------

validcell => soilvars(grid)%validcell

bulk   => soilvars(grid)%bulk
rock   => soilvars(grid)%rock
OrgM   => soilvars(grid)%OrgM
Vsand  => soilvars(grid)%Vsand
Vclay  => soilvars(grid)%Vclay
Vsilt  => soilvars(grid)%Vsilt
VOrgM  => soilvars(grid)%VOrgM
Csolid => soilvars(grid)%Csolid
Ksolid => soilvars(grid)%Ksolid
Kdry   => soilvars(grid)%Kdry

whc    => soilvars(grid)%whc
Ksat   => soilvars(grid)%Ksat
Tsat   => soilvars(grid)%Tsat
Tfield => soilvars(grid)%Tfield
Twilt  => soilvars(grid)%Twilt
Psat   => soilvars(grid)%Psat
Bexp   => soilvars(grid)%Bexp
Wliq   => soilvars(grid)%Wliq
Wice   => soilvars(grid)%Wice
Tsoil  => soilvars(grid)%Tsoil
Tsoiln => soilvars(grid)%Tsoiln
Tliq   => soilvars(grid)%Tliq
Tice   => soilvars(grid)%Tice
Psi    => soilvars(grid)%Psi

!-------------------------

allocate(valid(nl))

valid = .true.

! Assume 1% of soil organic matter (missing data from soil grid?)
OrgM = 1.

! Assume typical soiltype
soiltype = 1

bulkloop : do l = 1, nl

  sand  = soilvars(grid)%sand(l)
  clay  = soilvars(grid)%clay(l)
  cfvo  = soilvars(grid)%cfvo(l)
  OM    = soilvars(grid)%OrgM(l)

  silt = 100. - (sand + clay)

  if (sand < 0. .or. sand == missing_sp) then
    valid(l) = .false.
    cycle
  end if

  ! Conversion factor from organic matter to organic carbon
  Csoil = OrgM(l) / omcf  !layer

  orgC = Csoil

  ! Because bulk density depends strongly on organic matter content and
  ! weakly on wilting point water content, we guess an initial value and
  ! iterate to a stable solution
  it = 1                 ! Iteration number

  T1500 = 0.1

  bulk(l) = fbulk(orgC,clay,silt,T1500*100.,zpos(l))

  !-----

  do

    ! Convert SOM to mass fraction and calculate the difference in OM content
    soilmass = bulk(l) * 1.e6 * dz(l)                    !(g) (m-2)  note: g cm-3 * (m*100=cm) 100cm * 100cm = g

    orgC = max(100. * Csoil / soilmass,0.)            !(mass %)

    ! Recalculate bulk
    blk0 = fbulk(orgC,clay,silt,T1500*100.,zpos(l))   !units (g cm-3)

    ! Calculate wilting point, field capacity, and saturation, needs input in fractions not percent
    OM = orgC * omcf                  !(mass %)

    if (OM >= 30.) soiltype = 3       !humic soil

    call calctheta(sand/100.,clay/100.,OM/100.,blk0,Tsat(l),T33,T1500,soiltype)

    ! Recalculate bulk
    bulk(l) = fbulk(orgC,clay,silt,T1500*100.,zpos(l))  !units (g cm-3)

    if (abs(bulk(l) - blk0) < 0.001 .or. it > 50) exit

    ! Save bulk density value from latest iteration
    blk0 = bulk(l)

    it = it + 1

  end do

  !------

  ! With the final value for bulk density, recalculate porosity
  call calctheta(sand/100.,clay/100.,OM/100.,bulk(l),Tsat(l),T33,T1500)

  !------

  ! Stored bulk density for soil layer
  bulk(l) = bulk(l) * 1000.                       ! Convert g cm-3 to kg m-3

  ! Calculate course fragment content by mass
  rock(l) = fRock(cfvo/100., bulk(l))             ! mass by fraction

  ! Save soil water volumetric water content at saturation
  Tsat(l) = Tsat(l)                               ! fraction / m3 m-3

  ! Save soil water voulumetric water content at field capacity and wilting point
  Tfield(l) = T33                                 ! fraction / m3 m-3
  Twilt(l)  = T1500                               ! fraction / m3 m-3

  ! Update layer-integrated WHC = Tfield - Twilt
  whc(l) = 1000. * dz(l) * (T33 - T1500)           ! mm

  ! Calculate hydraulic condictivity at saturation
  Ksat(l) = fKsat(Tsat(l),T33,T1500)
  Ksat(l) = Ksat(l) / 3600.                       ! Convert mm h-1 to mm s-1

  ! Calculate soil water matric potential at saturation
  Psat(l) = fPsat(sand*(1.-rock(l)))              ! mm

  ! Calculate B exponent used in the Brooks & Corey Pedotransfer Functions
  ! NOTE: Reduce the sand,silt,clay mass percent by the coarse fragment mass percent (while retaining the relative weight percents
  Bexp(l) = fBexp(clay*(1.-rock(l)), silt*(1.-rock(l)), sand*(1.-rock(l)), rock(l), OrgM(l))

  ! if (rock(l) /= 0.) print *, it, Tsat(l), Tfield(l), whc(l)

end do bulkloop

! if (rock(1) /= 0.) print *, sum(Tsat(1:5)) / 5., sum(Tfield(1:5)) / 5., sum(Twilt(1:5)) / 5., sum(whc(1:5)) / 1000.

!-------------------------

! Initial soil temperature preset and calculate dry soil thermal conductivity and heat capacity
conductivityloop : do l = 1, nl

  !inital presets
  Tsoil(l) = Tfreeze + 5. !set to 5 deg. C.
  Tsoiln(l) = Tfreeze + 5.05

  !------

  OM   = 1.
  sand = soilvars(grid)%sand(l)
  clay = soilvars(grid)%clay(l)
  cfvo = soilvars(grid)%cfvo(l) * 0.01      ! Convert from percent to fraction
  silt = 100. - (sand + clay)

  !find the volume fractions of the soil consituents.
  VOrgM(l) = OM   * 0.01 / omd
  Vsand(l) = sand * 0.01 / sandd
  Vclay(l) = clay * 0.01 / clayd
  Vsilt(l) = silt * 0.01 / siltd

  Vtotal = VOrgM(l) + Vsand(l) + Vclay(l) + Vsilt(l) + cfvo  ! Total volume including course fragments

  VOrgM(l) = VOrgM(l) / Vtotal
  Vsand(l) = Vsand(l) / Vtotal
  Vclay(l) = Vclay(l) / Vtotal
  Vsilt(l) = Vsilt(l) / Vtotal

  ! if (sand > 0.) print *, VOrgM(l), Vsand(l), Vclay(l), Vsilt(l), cfvo, Vtotal

  !------

  ! Calculate soil solids volumetric heat capcity at layer midpoint (J m-3 K-1)
  Csolid(l) = VOrgM(l) * Corg + Vsand(l) * Csand + Vclay(l) * Cclay + cfvo * Crock + Vsilt(l) * Csilt

    Vquartz = sandquartz * Vsand(l) + rockquartz * cfvo

    interm  = Kquartz**Vquartz * Kmineral**(1. - VOrgM(l) - Vquartz)

  ! Calculate soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
  Ksolid(l) = Kom**VOrgM(l) * interm ! B&A Eqn 15

    Kdrysolid(l) = Kom_dry**VOrgM(l) * interm ! W m-1 K-1; B&A Eqn 15

  ! if (l <= gnl+1) then  !for soil column

    densp = 1. / (VOrgM(l) + (1. - OM * 0.01) / clayd) ! Eqn 3

  ! Calculate dry soil thermal conductivity at layer midpoint (W m-1 K-1)
  Kdry(l) = ((a * Kdrysolid(l) - Kair) * bulk(l) * 0.001 + Kair * densp) &
            / (densp - (1. - a) * bulk(l) * 0.001)   ! Eqn 16

  ! else  !set bedrock thermal conductivity
  !
  !   Kdry(l) = bedrock_tc
  !
  ! end if

  ! if (sand > 0.) print *, Ksolid(l), Csolid(l), Kdry(l), bulk(l)

end do conductivityloop

!-------------------------

! Initialize hydrology initial pre-sets for first run.
hydroloop : do l = 1, nl

  Tliq(l) = 0.3 * Tsat(l)

  ! peatlands get saturated soil columns
  ! if (peat) Tliq(l) = Tsat(l)  !FLAG test.  ---> peat dsiabled (Leo Lai Jul 2021)

  Wliq(l) = Tliq(l) * dz(l) * pliq

  ! Assume no ice.
  Tice(l) = 0.
  Wice(l) = 0.

  Psi(l) = Psat(l) * (Tliq(l) / Tsat(l))**(-Bexp(l))

  if (Psi(l) < -1.e8) Psi(l) = -1.e8


end do hydroloop

!-------------------------

! assign a typical rock value for non-soil layers
! whc  = 0.03
! Ksat = 1.e-3
! valid = .true.

validcell = .true.

if (.not. valid(1)) validcell = .false.

if (.not. validcell) then

  whc  = missing_sp
  Ksat = missing_sp
  Tsat = missing_sp
  Psat = missing_sp
  Tfield = missing_sp
  Twilt = missing_sp
  Bexp = missing_sp

  Vsand = missing_sp
  Vclay = missing_sp
  Vsilt = missing_sp
  VOrgM = missing_sp
  rock = missing_sp
  bulk = missing_sp

  Csolid = missing_sp
  Ksolid = missing_sp
  Kdry = missing_sp

  Wliq = missing_sp
  Wice = missing_sp
  Tsoil = missing_sp
  Tsoiln = missing_sp
  Tliq = missing_sp
  Tice = missing_sp
  Psi = missing_sp

end if

!-------------------------

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

! if (sum(whc(1:5)) > 0.) print *, sum(whc(1:5)) / 10.
! if (soilvars(grid)%whc(1) /= missing_sp) print *, Tsat(1:5)


end subroutine soilprep

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

function fbulk(orgC,clay,silt,W15,depth)

!equations from Heuscher et al., (2005) Soil Sci. Soc. Am. J. 69

real(sp) :: fbulk  !bulk density (g cm-3)

real(sp), intent(in) :: orgC   !mass %
real(sp), intent(in) :: clay   !mass %
real(sp), intent(in) :: silt   !mass %
real(sp), intent(in) :: W15    !water content at -15 bar (-1500 kPa) (mass %)
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

function fRock(Vcf,bulk)

!correction factor for rock fragments (from Skirvin correction to WEPP model)
!this function takes the volume fraction and converts it to a mass fraction for
!coarse fragments.

! Function copied from ARVE-DGVM original code (Leo O Lai, Jul 2021)

implicit none

!real(dp) :: fVcf
!real(dp) :: rock
real(sp) :: Vcf
real(sp) :: frock
real(sp) :: bulk

real(sp), parameter :: a = 2685.
!---

  if (Vcf > 0.) then
    frock = a / (bulk / Vcf - bulk + a)
  else
    frock = 0.
  end if

  !depricated function., presently not required. JM 01.06.2010.
  !fVcf = (rock * 0.01 * bulk) / (rock * 0.01 * bulk + a * (1. - rock * 0.01))
  !Alternative method to estimate from Brakensiek & Rawls 1994 Catena
  !this method is not presently used.
  !frock = 2. * Vcf / (1. + Vcf)

end function fRock

!---------------------------------------------------------------------------

function fKsat(Tsat,T33,T1500)

!equations from Saxton & Rawls (2006) Soil Sci. Soc. Am. J. 70

real(sp) :: fKsat  !(mm h-1)

real(sp), intent(in) :: Tsat      ! fraction
real(sp), intent(in) :: T33       ! fraction
real(sp), intent(in) :: T1500     ! fraction

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

function fPsat(sand)

! Function to calculate saturated soil matric potential
! Closely and directly related to largest pore size
! Relation from CLM4 eqn 7.87

! Function copied from ARVE-DGVM original code (Leo O Lai, Jul 2021)

implicit none

real(sp) :: fPsat
real(sp), intent(in) :: sand   ! Soil sand content by mass (percent)

real(sp), parameter :: a = 1.8800
real(sp), parameter :: b = 0.0131

fPsat = -10. * 10. ** (a - b * sand)

end function fPsat

!---------------------------------------------------------------------------

function fBexp(clay,silt,sand,rock,sorg)

! Based upon Bloemen, 1980 Z. Pflanzenernaehr. Bodenkd. 143, 581-605
! Since Bloemen assumes a lower limit of 1.4 that value is used to
! convert the 'n' value to a Bexp value suitable for use in the rest of the model
! CLM uses a lower limit of '2'.

! Function copied from ARVE-DGVM original code (Leo O Lai, Jul 2021)

implicit none

real(sp) :: fBexp
real(sp), intent(in) :: clay   ! Soil clay content (percent)
real(sp), intent(in) :: silt   ! Soil clay content (percent)
real(sp), intent(in) :: sand   ! Soil sand content (percent)
real(sp), intent(in) :: rock   ! Soil rock content (mass fraction)
real(sp), intent(in) :: sorg   ! Soil organic content (percent)

real(sp), dimension(4), parameter :: sclass = [ 1., 26., 1125., 126000. ]  !size classes in um 0-2,2-50,50-2000,2000-250000 (~6000) so took means
real(sp), parameter :: b = 1.4
real(sp), parameter :: c = 4.536
real(sp), parameter :: d = 0.75
real(sp), parameter :: e = 1.6

real(sp), dimension(5) :: wps
real(sp), dimension(4) :: wps_cum
real(sp), dimension(3) :: tg
real(sp), dimension(4) :: f
real(sp) :: fsum
real(sp) :: n
integer :: i
integer :: base

!---

wps(1) = clay
wps(2) = silt
wps(3) = sand
wps(4) = rock * 100.
wps(5) = sorg

!NOTE: if there is no clay in the soil, then you can get a circumstance of divide by zero in equation 8
! to prevent that, tell it to simply skip the clay (moves the start of the loops to the second consitituent, which is silt.
! JM 23.04.2010

 if (wps(1) > 0.) then
   base = 1
 else
   base = 2
 end if


 do i = base,4 !make the Pi in the eqn 8

  wps_cum(i) = sum(wps(1:i))

 end do

 do i = base,3 !Makes the tgi (eqn 8) for the eqn 9

  tg(i) = log10(wps_cum(i+1) / wps_cum(i)) / log10(sclass(i+1) / sclass(i))

 end do

 do i = base,3

  f(i) = wps(i+1) * tg(i)

 end do

fsum = sum(f(1:3)) / sum(wps(2:4))

if (sorg > 0.) then  ! equation 12

  n = (b + c * (exp(fsum) - 1.))- d * (fsum**e) * log10(wps(5))
  fBexp = 3. / (n - 1.4)

else

  n = b + c * (exp(0.3 * fsum) - 1.)  !equation 11
  fBexp = 3. / (n - 1.4)

end if

end function fBexp

!---------------------------------------------------------------------------

end module soilstatemod
