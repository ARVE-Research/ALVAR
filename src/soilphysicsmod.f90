module soilphysicsmod

! Module adapted from ARVE-DGVM to calculate soil physics (Leo Lai, Aug 2021)

use parametersmod, only : i2,i4,sp,dp,missing_sp,Tfreeze

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine soilthermalprop(grid)

! This has been updated based upon Ballard & Arp 2005, J. Environ Eng. Sci 4:549-588 to
! account for organic matter and coarse fragments in the soil
! Code adapted from ARVE-DGVM (Leo Lai, Aug 2021) --> original in daily_calcs day/night loop

use metvarsmod,   only : soilvars
use soilstatemod, only : nl,dz,dzmm,zpos,zposmm,zipos,ziposmm

implicit none

integer(i4), intent(in) :: grid

! Pointers
logical , pointer :: validcell
real(sp), pointer, dimension(:) :: Tsat    ! Soil water volumetric water content at saturation (fraction)
real(sp), pointer, dimension(:) :: Bexp    ! Soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(sp), pointer, dimension(:) :: Psat    ! Soil water matric potential at saturation (mm)
real(sp), pointer, dimension(:) :: Tsoil   ! Soil temperature (K)
real(sp), pointer, dimension(:) :: Tliq    ! Soil liquid water content (fraction)
real(sp), pointer, dimension(:) :: Tice    ! Soil ice content (fraction)
real(sp), pointer, dimension(:) :: Vcf     ! FRACTION of coarse fragments (rocks) by volume
real(sp), pointer, dimension(:) :: Vom     ! Volumetric fraction of soil organic matter
real(sp), pointer, dimension(:) :: Vsand   ! Volumetric fraction of sand
real(sp), pointer, dimension(:) :: Kdry    ! Soil thermal conductivity of dry natural soil (W m-1 K-1)
real(sp), pointer, dimension(:) :: Ksolid  ! Soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
real(sp), pointer, dimension(:) :: Csolid  ! Soil solids volumetric heat capacity at layer midpoint (J m-3 K-1)
real(sp), pointer, dimension(:) :: Kl      ! Soil thermal conductivity across layer boundary (W m-1 K-1)
real(sp), pointer, dimension(:) :: Cl      ! Instantaneous soil volumetric heat capacity at layer midpoint (J m-3 K-1)
real(sp), pointer, dimension(:) :: Kh      ! Snow thermal conductivity at layer midpoint (W m-1 K-1)

! Parameters
real(sp), parameter :: alpha = 0.24_sp     ! +/- 0.04 constant for Knum calc.from Balland and Arp
real(sp), parameter :: beta  = 18.3_sp     ! +/- 1.1 constant for Knum calc.from Balland and Arp
real(sp), parameter :: Cliq  = 4.18800e3   ! Heat capacity of liquid water  (J kg-1 K-1)
real(sp), parameter :: Cice  = 2.11727e3   ! Heat capacity of ice (typical) (J kg-1 K-1)
real(sp), parameter :: Kliq  = 0.57_sp     ! Thermal conductivity of liquid water (typical) (W m-1 K-1)
real(sp), parameter :: Kice  = 2.2_sp      ! Thermal conductivity of ice (typical) (W m-1 K-1)
real(sp), parameter :: Kair  = 0.0243_sp   ! Thermal conductivity of dry air (W m-1 K-1)

! Local variables
real(sp) :: Knum          ! Kersten number (CLM technical note eqn 6.63, pg. 94)
real(sp) :: Ktsat         ! Soil saturated thermal conductivity (W m-1 K-1)
real(sp) :: Sr            ! Soil wetness with respect to saturation
real(sp) :: Knum1         ! Interm variables.
real(sp) :: Knum2
integer :: l

!-------------------------

validcell => soilvars(grid)%validcell

Tsat   => soilvars(grid)%Tsat
Bexp   => soilvars(grid)%Bexp
Psat   => soilvars(grid)%Psat
Tsoil  => soilvars(grid)%Tsoil
Tliq   => soilvars(grid)%Tliq
Tice   => soilvars(grid)%Tice
Vcf    => soilvars(grid)%cfvo
Vom    => soilvars(grid)%VOrgM
Vsand  => soilvars(grid)%Vsand
Kdry   => soilvars(grid)%Kdry
Ksolid => soilvars(grid)%Ksolid
Csolid => soilvars(grid)%Csolid
Kl     => soilvars(grid)%Kl
Cl     => soilvars(grid)%Cl
Kh     => soilvars(grid)%Kh

!-------------------------

! While the bottom layers (gnl+1,nl) are not hydrologically active, they do freeze/thaw
! so soil thermal properties are calculated for all soil levels.

! Thermal conductivity and heat capacity of soil solids (from Hillel, Env. Soil Phys Book)
do l = 1,nl

  ! Saturated thermal conductivity
  if (Tsoil(l) > Tfreeze) then !non-frozen soils

    Ktsat = Ksolid(l)**(1.0 - Tsat(l)) * Kliq**Tsat(l) !eqn 12 in Balland and Arp

  else  !frozen soils

    Ktsat = Ksolid(l)**(1.0 - Tsat(l)) * Kliq**Tliq(l) * Kice**(Tsat(l) - Tliq(l)) ! eqn. 13 in Balland & Arp

  end if

  !------
  ! Wetness with respect to saturation
  Sr = min((Tliq(l) + Tice(l)) / Tsat(l), 1.0)  !eqn. 6

  ! Kersten number
  if (Tsoil(l) > Tfreeze) then    !thawed soils (eqn 17 balland & arp)

    Knum1 = 0.5 * (1.0 + Vom(l) - alpha * Vsand(l) - Vcf(l) / 100.)

    Knum2 = ((1.0 / (1.0 + exp(-beta * Sr)))**3 - ((1.0 - Sr) / 2.0)**3)**(1.0 - Vom(l))

    Knum = Sr**Knum1 * Knum2

  else  ! Frozen or partially frozen soils   (eqn 18 balland & arp)

    Knum = Sr**(1.0 + Vom(l))

  end if

  !------
  ! Instantaneous soil thermal conductivity at layer midpoint
  Kh(l) = Knum * (Ktsat - Kdry(l)) + Kdry(l)    ! (eqn 5 balland & arp)

  ! Soil heat capacity at layer midpoint
  ! NOTE: new formulation - Apr 18 08 JM
  Cl(l) = ((Csolid(l) * (1.0 - Tsat(l))) + (Tice(l) * 1000. * Cice) + (Tliq(l) * 1000. * Cliq)) !+ &
  ! (1000._dp * Cair * (Tsat(l) - Tliq(l) - Tice(l))) ) !NOTE the last part is to account for the heat capacity of the air
  ! within the pore space, ignored at present.

end do

!-------------------------

do l = 1,nl

  ! Find soil thermal conductivity across layer boundary (W m-1 K-1)
  if (l < nl) then

    Kl(l) = (Kh(l) * Kh(l+1) * (zpos(l+1) - zpos(l))) /  &
            (Kh(l) * (zpos(l+1) - zipos(l)) + Kh(l+1) * (zipos(l) - zpos(l)))

  else

    Kl(l) = 0.0  ! No flux bottom condition

  end if

end do

! print *, Tsat(1), Cl(1), Kdry(nl), Ktsat, Kl


end subroutine soilthermalprop

!---------------------------------------------------------------------

subroutine resistance(grid,day)

! Calculation of resistance to sensible and latent heat transfer for
! vegetated and non-vegetated surfaces. Based upon the CLM 3.1 formulation with multiple modifications
! NOTE: there is presently no difference in the calculation of sensible vs. latent resistances!
! Added  to ARVE-DGVM Oct 08 08 Joe Melton 2008, adapted for ALVAR (Leo Lai Aug 2021)
!     --> Currently simplied version with ONLY BARE SURFACE

use metvarsmod,   only : dayvars,soilvars,gprint,lprint

integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

! Pointer variables
real(sp), pointer :: wind                      ! Wind speed (m s-1)
real(sp), pointer :: raw_b                     ! Aerodynamic resistance for latent heat transfer btwn bare ground & atmosphere (s m-1)
real(sp), pointer :: rah_b                     ! Aerodynamic resistance for sensible heat transfer btwn bare ground & atmosphere (s m-1)

! Parameters
real(sp), parameter :: refH = 5._sp            ! reference height (m)
real(sp), parameter :: vkarm = 0.4_sp          ! von Karman constant (unitless)
real(sp), parameter :: visco = 1.5e-5          ! kinematic viscosity of air (m2 s-1)
real(sp), parameter :: z0mg  = 0.01_sp         ! momentum roughness length for bare soil (m, CLM eqn 3.49)
real(sp), parameter :: a = 0.13_dp

! Local variables
real(sp) :: mroughl        ! Momentum calculation roughness length (cm)
real(sp) :: roughl_b       ! Roughness length of bare surface
real(sp) :: ustar_b        ! Friction velocity estimated from scalar wind speed (m s-1)
real(sp) :: z_b            ! Height above bare surface (cm)
real(sp) :: zpdisp_b       ! Zero plane dispalcement of bare surface (cm)

!-------------------------

wind  => dayvars(grid,day)%wind

raw_b => soilvars(grid)%raw_b
rah_b => soilvars(grid)%rah_b

!-------------------------

mroughl = z0mg

! Find the resistances for bare ground

z_b = refH          ! 5m in m
zpdisp_b = 0._sp    ! Zero plane displacement (m)

ustar_b = 0.14 * wind   ! Friction velocity estimated from the scalar wind speed (Weber 1999 Boundary Layer Meteor, 93, 197)

! Convert from momentum roughness length to sensible and latent roughness length. This is because the transfer
! of momentum is affected by pressure fluctuations in the turbulent waves behind the roughness elements, while for
! heat and water vapour transfer no such dynamical mechanism exists.(m)

roughl_b = mroughl * exp(-a * (ustar_b * mroughl / visco)**0.45)  !CLM 5.67

! The aerodynamic resistance for momentum transfer between the ground and the atmosphere (s m-1)
! from Szeicz et al. 1969 Water Resources Res. (5) 380. No difference is accounted for between
! sensible and latent heat transfers
rah_b = (log((z_b - zpdisp_b) / roughl_b))**2 / (vkarm**2 * wind) !s m-1
raw_b = rah_b

if (wind == 0) rah_b = 4000.
if (wind == 0) raw_b = 4000.

! if (grid==gprint .and. lprint) print *, wind, rah_b, raw_b

end subroutine resistance

!---------------------------------------------------------------------

subroutine soiltemperature(grid,day,i)

use metvarsmod,   only : dayvars,soilvars,gprint,lprint
use soilstatemod, only : nl,dz,dzmm,zpos,zposmm,zipos,ziposmm
use utilitiesmod, only : tridiag

integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day
integer(i4), intent(in) :: i

! Pointer variables
integer(i4), pointer :: dayhour       ! day time hours (h) --> from current day sunrise to sunset
integer(i4), pointer :: nighthour     ! night time hours (h) --> from current day sunset to next day sunrise
real(sp), pointer :: Patm             ! Atmospheric pressure to sea-level (Pa)
real(sp), pointer :: tday             ! Mean daytime temperature (degC)
real(sp), pointer :: tnight           ! Mean nighttime temperature (degC)
real(sp), pointer :: srad             ! downwelling surface shortwave radiation (kJ m-2 d-1)
real(sp), pointer :: srad_dir         ! direct beam downwelling shortwave raditaion (kJ m-2 d-1)
real(sp), pointer :: srad_dif         ! diffuse downwelling shortwave raditaion (kJ m-2 d-1)
real(sp), pointer :: lrad             ! upswelling surface longwave radiation (kJ m-2 d-1)
real(sp), pointer :: dayl             ! Daylength (h)
real(sp), pointer, dimension(:) :: Tsoil       ! Soil temperature (K)
real(sp), pointer, dimension(:) :: Tsoiln      ! Soil temperature for precious timestep (K)
real(sp), pointer, dimension(:) :: Kl          ! Soil thermal conductivity across layer boundary (W m-1 K-1)
real(sp), pointer, dimension(:) :: Cl          ! Instantaneous soil volumetric heat capacity at layer midpoint (J m-3 K-1)
real(sp), pointer, dimension(:) :: Kh          ! Snow thermal conductivity at layer midpoint (W m-1 K-1)
real(sp), pointer :: raw_b
real(sp), pointer :: rah_b

real(sp), parameter :: Cair  = 1.25e3        ! Heat capacity of dry air (J kg-1 K-1) (CLM parameterization)
real(sp), parameter :: Rspec = 287.058       ! Specific gas constant of dry air
real(sp), parameter :: Tstune = 0.34_sp      ! Tuning factor to turn first layer T into surface T (CLM parameterization, pg 88).
real(sp), parameter :: CNfac = 0.5_sp        ! Crank Nicholson factor between 0 and 1

! Local variables
real(sp) :: Swg_b       ! Shortwave radiation absorbed by bare ground (W m-2)
real(sp) :: Lwg_b       ! Longwave radiation absorbed by bare ground (W m-2)
! real(sp) :: rah_b       ! Resistance to sensible heat transfer (bare ground and atm)(s m-1)
! real(sp) :: raw_b       ! Resistance to latent heat transfer (bare ground and atm)(s m-1)
real(sp) :: Hg_b        ! Sensible heat flux into the surface (W m-2) bare surfaces
real(sp) :: Eg_b        ! Latent heat flux into the surface (mm s-1)(bare surfaces)
real(sp) :: dHgdt_b     ! Derivative of the sensible heat flux into the surface (W m-2) w.r.t temp (bare surf)
real(sp) :: dEgdt_b     ! Derivative of the latent heat flux into the surface (W m-2) w.r.t temp (bare surf)
real(sp) :: Datm        ! Atmospheric density (kg m-3)
real(sp) :: hs          ! Net energy flux into the surface (W m-2)

! Tridiagonal solver variables
real(sp) :: delZtop
real(sp) :: dhsdT
real(sp) :: dzp
real(sp) :: dzm
real(sp), dimension(nl) :: fact
real(sp), dimension(nl) :: Fhti

real(dp), dimension(nl) :: avect             !vectors for tridiagonal solver
real(dp), dimension(nl) :: bvect
real(dp), dimension(nl) :: cvect
real(dp), dimension(nl) :: rvect
real(dp), dimension(nl) :: Tnsoi

integer :: dt           ! Time-step length (seconds)
integer :: a
integer :: b
integer :: l
integer :: h
integer :: hr

!-------------------------

dayhour   => dayvars(grid,day)%dayhour
nighthour => dayvars(grid,day)%nighthour
Patm      => dayvars(grid,day)%Patm
tday      => dayvars(grid,day)%tday
tnight    => dayvars(grid,day)%tnight
srad      => dayvars(grid,day)%srad
srad_dir  => dayvars(grid,day)%srad_dir
srad_dif  => dayvars(grid,day)%srad_dif
lrad      => dayvars(grid,day)%lrad
dayl      => dayvars(grid,day)%dayl

Tsoil     => soilvars(grid)%Tsoil
Tsoiln    => soilvars(grid)%Tsoiln
Kl        => soilvars(grid)%Kl
Cl        => soilvars(grid)%Cl
Kh        => soilvars(grid)%Kh
raw_b     => soilvars(grid)%raw_b
rah_b     => soilvars(grid)%rah_b

!-------------------------

Eg_b = 0.0
dEgdt_b = 0.0

! Calculate sensible heat flux between the air and soil surface (CLM eqn 5.60, CLM 4.0 Eqn 5.61)
! NOTE: Potential temperature of atmosphere substituted for simple surface air temperature

if (i == 1) then          ! Day time

  Datm = Patm / (Rspec * (tday + Tfreeze))

  Hg_b = -Datm * Cair * (tday + Tfreeze - Tsoil(1)) / rah_b
  dHgdt_b = Datm * Cair / rah_b

  Swg_b = srad_dir * 1000. / (dayl * 3600.)
  Lwg_b = lrad * 1000. / (dayl * 3600.)

else if (i == 2) then     ! Night time

  Datm = Patm / (Rspec * (tnight + Tfreeze))

  Hg_b = -Datm * Cair * (tnight + Tfreeze - Tsoil(1)) / rah_b
  dHgdt_b = Datm * Cair / rah_b

  Swg_b = 0.0
  Lwg_b = 0.0

end if

hs = Swg_b - Lwg_b - Hg_b

dhsdT = dHgdt_b

!-------------------------
! Compile met variables based on dayhour / nighthour

if (i == 1) then            ! Day time
  hr = dayhour
else if (i == 2) then       ! Night time
  hr = nighthour
end if

!-------------------------
! Time-step dt = 3600 for now as subroutine written for constant 1 hour time-step
dt = 3600.

!-------------------------
! Run water model in hourly time-step based on length of dayhour / nighthour
hourloop : do h = 1, hr

  ! Set-up tridiagonal system for current time-step soil temperature

  ! Top snow/soil layer is the layer-averaged temp, so to be more accurate
  ! we adjust the heat capacity of the top layer. This will give a more realistic top soil layer
  ! temperature. See P.107 of CLM 4.0. Eqn 6.29
  ! Start with top soil layer
  l = 1

  delZtop = 0.5_dp * (zpos(l) - zipos(l-1) + Tstune * (zpos(l+1) - zipos(l-1)))

  fact(l) = dt / Cl(l) * dz(l) / delZtop
  Fhti(l) = Kl(l) * (Tsoil(l+1) - Tsoil(l)) / (zpos(l+1) - zpos(l))

  ! Middle layers

  do l = 2, nl-1
    fact(l) = dt / Cl(l)
    Fhti(l) = -Kl(l) * (Tsoil(l) - Tsoil(l+1)) / (zpos(l+1) - zpos(l))
  end do

  ! Bottom soil layer

  fact(nl) = dt / Cl(nl)
  Fhti(nl) = 0._dp

  ! Set up the coefficients for the tridiagonal matrix solver
  ! top snow/soil layer
  l = 1

    dzp = zpos(l+1) - zpos(l)

    avect(l) = 0._sp
    bvect(l) = 1._sp + (1._sp - CNfac) * fact(l) * (Kl(l) / dzp) - fact(l) * dhsdT
    cvect(l) =       - (1._sp - CNfac) * fact(l) *  Kl(l) / dzp
    rvect(l) = Tsoil(l) + fact(l) * (hs - dhsdT * Tsoil(l) + CNfac * Fhti(l))

  ! Middle soil layers
  do l = 2, nl-1

    dzm = zpos(l) - zpos(l-1)
    dzp = zpos(l+1) - zpos(l)

    avect(l) =       - (1._dp - CNfac) * fact(l) *  Kl(l-1) / dzm
    bvect(l) = 1._dp + (1._dp - CNfac) * fact(l) * (Kl(l)   / dzp + Kl(l-1) / dzm)
    cvect(l) =       - (1._dp - CNfac) * fact(l) *  Kl(l)   / dzp
    rvect(l) = Tsoil(l) + CNfac * fact(l) * (Fhti(l) - Fhti(l-1))

  end do

  ! bottom soil layer
  l = nl

    dzm = zpos(l) - zpos(l-1)

    avect(l) =       - (1._sp - CNfac) * fact(l) *  Kl(l-1) / dzm
    bvect(l) = 1._sp + (1._sp - CNfac) * fact(l) *  Kl(l-1) / dzm
    cvect(l) = 0._sp
    rvect(l) = Tsoil(l) - CNfac * fact(l) * Fhti(l-1)

  ! Call tridiagonal solver
  a = 1
  b = nl

  call tridiag(avect(a:b),bvect(a:b),cvect(a:b),rvect(a:b),Tnsoi(a:b))

  Tsoiln = Tsoil
  Tsoil = Tnsoi

end do hourloop

! if (lprint .and. grid==gprint) print *, day, i, hs, rah_b, tday,tnight, Tsoil-Tfreeze


end subroutine soiltemperature


end module soilphysicsmod
