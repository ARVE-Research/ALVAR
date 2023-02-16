module soilphysicsmod

! Module adapted from ARVE-DGVM to calculate soil physics (Leo Lai, Aug 2021)

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine soilthermalprop(validcell,snl,dz,zpos,zipos,Vcf,Vsand,Vom,Tsat,Psat,Bexp,Tsoil,Tliq,Tice,Wliq,Wice,Wsno,&
                           Csolid,Ksolid,Kdry,Kl,Cl,Kh)

! This has been updated based upon Ballard & Arp 2005, J. Environ Eng. Sci 4:549-588 to
! account for organic matter and coarse fragments in the soil
! Code adapted from ARVE-DGVM (Leo Lai, Aug 2021) --> original in daily_calcs day/night loop

use parametersmod, only : i4,sp,Tfreeze
use statevarsmod,  only : nl,ns

implicit none

logical,                intent(in)    :: validcell
integer(i4),            intent(in)    :: snl
real(sp), dimension(ns:nl), intent(in)    :: dz
real(sp), dimension(ns:nl), intent(in)    :: zpos         ! Midpoint z position (depth) of soil layer (m)
real(sp), dimension(ns:nl+1), intent(in)    :: zipos        ! Snow/soil layer interface z position (m)
real(sp), dimension(:), intent(inout) :: Vcf          ! FRACTION of coarse fragments (rocks) by volume
real(sp), dimension(:), intent(inout) :: Vsand        ! Volumetric fraction of sand
real(sp), dimension(:), intent(inout) :: Vom          ! Volumetric fraction of soil organic matter
real(sp), dimension(ns:nl), intent(inout) :: Tsat         ! Soil water volumetric water content at saturation (fraction)
real(sp), dimension(ns:nl), intent(inout) :: Psat         ! Soil water matric potential at saturation (mm)
real(sp), dimension(ns:nl), intent(inout) :: Bexp         ! Soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(sp), dimension(ns:nl), intent(inout) :: Tsoil        ! Soil temperature (K)
real(sp), dimension(ns:nl), intent(inout) :: Tliq         ! Soil liquid water content (fraction)
real(sp), dimension(ns:nl), intent(inout) :: Tice         ! Soil ice content (fraction)
real(sp), dimension(ns:nl), intent(inout) :: Wliq           ! Soil liquid water content at layer midpoint (mm)
real(sp), dimension(ns:nl), intent(inout)  :: Wice       ! Soil ice content at layer midpoint (mm)
real(sp)                , intent(in) :: Wsno       ! Snow water equivalent of the snowpack (mm)
real(sp), dimension(nl), intent(inout) :: Csolid       ! Soil solids volumetric heat capacity at layer midpoint (J m-3 K-1)
real(sp), dimension(nl), intent(inout) :: Ksolid       ! Soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
real(sp), dimension(nl), intent(inout) :: Kdry         ! Soil thermal conductivity of dry natural soil (W m-1 K-1)
real(sp), dimension(ns:nl), intent(inout) :: Kl           ! Soil thermal conductivity across layer boundary (W m-1 K-1)
real(sp), dimension(ns:nl), intent(inout) :: Cl           ! Instantaneous soil volumetric heat capacity at layer midpoint (J m-3 K-1)
real(sp), dimension(ns:nl), intent(inout) :: Kh           ! Snow thermal conductivity at layer midpoint (W m-1 K-1)

! Parameters
real(sp), parameter :: alpha = 0.24_sp     ! +/- 0.04 constant for Knum calc.from Balland and Arp
real(sp), parameter :: beta  = 18.3_sp     ! +/- 1.1 constant for Knum calc.from Balland and Arp
real(sp), parameter :: Cliq  = 4.18800e3   ! Heat capacity of liquid water  (J kg-1 K-1)
real(sp), parameter :: Cice  = 2.11727e3   ! Heat capacity of ice (typical) (J kg-1 K-1)
real(sp), parameter :: pice  = 917.0       ! Density of ice (kg m-3)
real(sp), parameter :: Kliq  = 0.57_sp     ! Thermal conductivity of liquid water (typical) (W m-1 K-1)
real(sp), parameter :: Kice  = 2.2_sp      ! Thermal conductivity of ice (typical) (W m-1 K-1)
real(sp), parameter :: Kair  = 0.0243_sp   ! Thermal conductivity of dry air (W m-1 K-1)

! Local variables
integer(i4) :: l
real(sp) :: Knum          !Kersten number (CLM technical note eqn 6.63, pg. 94)
real(sp) :: Ktsat         !soil saturated thermal conductivity (W m-1 K-1)
real(sp) :: psno          !snow bulk density
real(sp) :: Sr            !soil wetness with respect to saturation
real(sp) :: Knum1         !interm variables.
real(sp) :: Knum2
real(sp) :: term_a        !Balland and Arp snow thermal conductivity adjustable term

!-------------------------

! While the bottom layers (gnl+1,nl) are not hydrologically active, they do freeze/thaw
! so soil thermal properties arescalculated for all soil levels.

!thermal conductivity and heat capacity of soil solids (from Hillel, Env. Soil Phys Book) and snow
do l = snl+1,nl

  !snow layers instantaneous conductivity and heat capacity
  if (l < 1) then

            !snow bulk density (kg m-3) eqn 6.66 CLM
            psno = (Wice(l) + Wliq(l)) / dz(l)

            ! snow thermal conductivity CLM original formulation
!            Kh(l) = Kair + (7.75e-5 * psno + 1.105e-6 * psno * psno) * (Kice - Kair) !(W m-1 K-1) eqn 6.65 CLM

            ! snow thermal conductivity from Balland and Arp 2005.
            term_a = 0.3
            Kh(l) = ((term_a * Kice - Kair) * psno * 1.e3 + Kair * pice * 1.e3) / (pice * 1.e3 - (1. - term_a) * psno * 1.e3)

            !snow volumetric heat capacity
            Cl(l) = (Wice(l) / dz(l)) * Cice + (Wliq(l) / dz(l)) * Cliq  !(J m-3 K-1) eqn 6.69 CLM

  else     !all soil layers

          !saturated thermal conductivity
          if (Tsoil(l) > Tfreeze) then !non-frozen soils
             Ktsat = Ksolid(l)**(1. - Tsat(l)) * Kliq**Tsat(l) !eqn 12 in Balland and Arp
          else  !frozen soils
             Ktsat = Ksolid(l)**(1. - Tsat(l)) * Kliq**Tliq(l) * Kice**(Tsat(l) - Tliq(l)) ! eqn. 13 in Balland & Arp
          end if

            !wetness with respect to saturation
            Sr = min((Tliq(l) + Tice(l)) / Tsat(l), 1.0)  !eqn. 6

            !Kersten number
            if (Tsoil(l) > Tfreeze) then    !thawed soils (eqn 17 balland & arp)

              Knum1 = (0.5 * (1. + Vom(l) - (alpha * Vsand(l)) - Vcf(l) / 100.))
              Knum2 = (((1. / (1. + exp(-beta * Sr)))**3.) - ((1. - Sr) / 2.)**3.)**(1. - Vom(l))

              Knum = Sr**Knum1 * Knum2

            else  !frozen or partially frozen soils   (eqn 18 balland & arp)

              Knum = Sr**(1. + Vom(l))

            end if

            !instantaneous soil thermal conductivity at layer midpoint
            Kh(l) = Knum * (Ktsat - Kdry(l)) + Kdry(l)    !(eqn 5 balland & arp)

            !soil heat capacity at layer midpoint
            !NOTE: new formulation - Apr 18 08 JM
            Cl(l) = ((Csolid(l) * (1. - Tsat(l))) + (Tice(l) * 1000. * Cice) + (Tliq(l) * 1000. * Cliq)) !+ &
                !(1000. * Cair * (Tsat(l) - Tliq(l) - Tice(l))) ) !NOTE the last part is to account for the heat capacity of the air
            !within the pore space, ignored at present.

  end if

end do

! Find soil/snow thermal conductivity across layer boundary (W m-1 K-1)
do l = snl+1,nl
  if (l < nl) then
    Kl(l) = (Kh(l) * Kh(l+1) * (zpos(l+1) - zpos(l))) /  &
             (Kh(l) * (zpos(l+1) - zipos(l+1)) + Kh(l+1) * (zipos(l+1) - zpos(l)))
  else
    Kl(l) = 0.  !no flux bottom condition
  end if
end do


!set the volumetric heat capacity of a snow layer if it is the inital one.
if (snl+1 == 1 .and. Wsno > 0. .and. dz(snl+1) /= 0.) then

  Cl(snl+1) = (Cl(snl+1) + Cice * Wsno / dz(snl+1))

  !FLAG not sure if these last 2 are necessary, I am adding them to cover my bases- JM 15.06.2010
  !I put this back to snl+1, it was l, from cvs version 1.9 though it does not appear to make sense as
  ! l. JM 25.10.10

  !snow bulk density (kg m-3) eqn 6.66 CLM
  psno = (Wice(snl+1) + Wliq(snl+1)) / dz(snl+1)

  !snow thermal conductivity
  Kh(snl+1) = Kair + (7.75e-5 * psno + 1.105e-6 * psno * psno) * (Kice - Kair) !(W m-1 K-1) eqn 6.65 CLM

end if

! print *, Kl(5), Kh(5:6), Kdry(6), Ksolid(5), Ktsat, Sr, Tsat(6), Knum1, Knum2, Knum, Vom,Vcf(6),Vom(6)


end subroutine soilthermalprop

!---------------------------------------------------------------------

subroutine soilresistance(wind,raw_b,rah_b)

! Calculation of resistance to sensible and latent heat transfer for
! vegetated and non-vegetated surfaces. Based upon the CLM 3.1 formulation with multiple modifications
! NOTE: there is presently no difference in the calculation of sensible vs. latent resistances!
! Added  to ARVE-DGVM Oct 08 08 Joe Melton 2008, adapted for ALVAR (Leo Lai Aug 2021)
!     --> Currently simplied version with ONLY BARE SURFACE

use parametersmod, only : sp

implicit none

real(sp), intent(in)    :: wind         ! Wind speed (m s-1)
real(sp), intent(inout) :: raw_b        ! Aerodynamic resistance for latent heat transfer btwn bare ground & atmosphere (s m-1)
real(sp), intent(inout) :: rah_b        ! Aerodynamic resistance for sensible heat transfer btwn bare ground & atmosphere (s m-1)

! Parameters
real(sp), parameter :: refH  = 5._sp            ! reference height (m)
real(sp), parameter :: vkarm = 0.4_sp          ! von Karman constant (unitless)
real(sp), parameter :: visco = 1.5e-5          ! kinematic viscosity of air (m2 s-1)
real(sp), parameter :: z0mg  = 0.01_sp         ! momentum roughness length for bare soil (m, CLM eqn 3.49)
real(sp), parameter :: a     = 0.13_sp

! Local variables
real(sp) :: mroughl        ! Momentum calculation roughness length (cm)
real(sp) :: roughl_b       ! Roughness length of bare surface
real(sp) :: ustar_b        ! Friction velocity estimated from scalar wind speed (m s-1)
real(sp) :: z_b            ! Height above bare surface (cm)
real(sp) :: zpdisp_b       ! Zero plane dispalcement of bare surface (cm)

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


end subroutine soilresistance

!---------------------------------------------------------------------

subroutine soiltemperature(i,snl,Patm,dayl,dayhour,nighthour,tday,tnight,srad_dir,lrad,dz,dzmm,zpos,zipos,&
                           Wliq,Wice,Tsat,Bexp,Psat,Tpor,Tliq,Tice,fice0,Wsno,zsno,qsnomelt,fsnow,Kl,Cl,Kh, &
                           fact,Fhti,ithaw,raw_b,rah_b,hs,dhsdT,Tsoil,Tsoiln)

use parametersmod, only : i4,sp,dp,Tfreeze
use statevarsmod,  only : nl,ns
use utilitiesmod,  only : tridiag
use soilphasechgmod, only : soilphasechg

implicit none

integer(i4),                     intent(in)    :: i
integer(i4),                     intent(in)    :: snl
real(sp),                        intent(in)    :: Patm          ! Atmospheric pressure to sea-level (Pa)
real(sp),                        intent(inout) :: dayl          ! Daylength (h)
integer(i4),                     intent(inout)    :: dayhour       ! Day time hours (h) --> from current day sunrise to sunset
integer(i4),                     intent(inout)    :: nighthour     ! Night time hours (h) --> from current day sunset to next day sunrise
real(sp),                        intent(in)    :: tday          ! Mean daytime temperature (degC)
real(sp),                        intent(in)    :: tnight        ! Mean nighttime temperature (degC)
real(sp),                        intent(in)    :: srad_dir      ! Direct beam downwelling shortwave raditaion (kJ m-2 d-1)
real(sp),                        intent(in)    :: lrad          ! Upswelling surface longwave radiation (kJ m-2 d-1)
real(sp),    dimension(ns:nl),   intent(in)    :: dz            ! Snow/soil layer thickness (m)
real(sp),    dimension(ns:nl),   intent(in)    :: dzmm          ! Thickness of the soil layers (mm)
real(sp),    dimension(ns:nl),   intent(in)    :: zpos          ! Midpoint z position (depth) of soil layer (m)
real(sp),    dimension(ns:nl+1), intent(in)    :: zipos         ! Snow/soil layer interface z position (m)
real(sp),    dimension(ns:nl),   intent(inout) :: Wliq          ! Soil liquid water content at layer midpoint (mm)
real(sp),    dimension(ns:nl),   intent(inout) :: Wice          ! Soil ice content at layer midpoint (mm)
real(sp),    dimension(ns:nl),   intent(inout) :: Tsat          ! Soil water volumetric water content at saturation (fraction)
real(sp),    dimension(ns:nl),   intent(inout) :: Bexp          ! Soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(sp),    dimension(ns:nl),   intent(inout) :: Psat          ! Soil water matric potential at saturation (mm)
real(sp),    dimension(ns:nl),   intent(inout) :: Tpor          ! Soil volumetric porosity (fraction)
real(sp),    dimension(ns:nl),   intent(inout) :: Tliq          ! Soil liquid water content (fraction)
real(sp),    dimension(ns:nl),   intent(inout) :: Tice          ! Soil ice content (fraction)
real(sp),    dimension(ns:nl),   intent(inout) :: fice0         ! Layer ice fraction, previous timestep
real(sp),                        intent(inout) :: Wsno          ! Snow water equivalent of the snowpack (mm)
real(sp),                        intent(inout) :: zsno          ! Total thickness of the snowpack (m)
real(sp),                        intent(inout) :: qsnomelt      ! Snow melt (kg m-2 s-1)
real(sp),                        intent(inout) :: fsnow         ! Fraction of the gridcell covered by snow (fraction)
real(sp),    dimension(ns:nl),   intent(inout) :: Kl            ! Soil thermal conductivity across layer boundary (W m-1 K-1)
real(sp),    dimension(ns:nl),   intent(inout) :: Cl            ! Instantaneous soil volumetric heat capacity at layer midpoint (J m-3 K-1)
real(sp),    dimension(ns:nl),   intent(inout) :: Kh            ! Snow thermal conductivity at layer midpoint (W m-1 K-1)
real(sp),    dimension(ns:nl),   intent(inout) :: fact          ! Factor used in computing tridiagonal coefficients
real(sp),    dimension(ns:nl),   intent(inout) :: Fhti          ! Heat flux across soil layer boundary (W m-2)
integer(i4), dimension(ns:nl),   intent(inout) :: ithaw
real(sp),                        intent(in)    :: raw_b         ! Aerodynamic resistance for latent heat transfer btwn bare ground & atmosphere (s m-1)
real(sp),                        intent(in)    :: rah_b         ! Aerodynamic resistance for sensible heat transfer btwn bare ground & atmosphere (s m-1)
real(sp),                        intent(inout)    :: hs            ! Net energy flux into the surface (W m-2)
real(sp),                        intent(inout)    :: dhsdT
real(sp),    dimension(ns:nl),   intent(inout) :: Tsoil         ! Soil temperature (K)
real(sp),    dimension(ns:nl),   intent(inout) :: Tsoiln        ! Soil temperature for previous timestep (K)

! Parameters
real(dp), parameter :: Cair  = 1.25e3        ! Heat capacity of dry air (J kg-1 K-1) (CLM parameterization)
real(dp), parameter :: Rspec = 287.058       ! Specific gas constant of dry air
real(dp), parameter :: Tstune = 0.34_dp      ! Tuning factor to turn first layer T into surface T (CLM parameterization, pg 88).
real(dp), parameter :: CNfac = 0.5_dp        ! Crank Nicholson factor between 0 and 1

! Local variables
real(sp) :: Swg_b       ! Shortwave radiation absorbed by bare ground (W m-2)
real(sp) :: Lwg_b       ! Longwave radiation absorbed by bare ground (W m-2)
real(sp) :: Hg_b        ! Sensible heat flux into the surface (W m-2) bare surfaces
real(sp) :: Eg_b        ! Latent heat flux into the surface (mm s-1)(bare surfaces)
real(sp) :: dHgdt_b     ! Derivative of the sensible heat flux into the surface (W m-2) w.r.t temp (bare surf)
real(sp) :: dEgdt_b     ! Derivative of the latent heat flux into the surface (W m-2) w.r.t temp (bare surf)
real(sp) :: Datm        ! Atmospheric density (kg m-3)
! real(sp) :: hs          ! Net energy flux into the surface (W m-2)

! Tridiagonal solver variables
real(sp) :: delZtop
! real(sp) :: dhsdT
real(sp) :: dzp
real(sp) :: dzm
real(dp), dimension(ns:nl) :: avect    ! Vectors for tridiagonal solver
real(dp), dimension(ns:nl) :: bvect
real(dp), dimension(ns:nl) :: cvect
real(dp), dimension(ns:nl) :: rvect
real(dp), dimension(ns:nl) :: Tnsoi

real(sp) :: dt           ! Time-step length (seconds)
real(sp) :: surf_atm_diff
integer(i4) :: a
integer(i4) :: b
integer(i4) :: l
integer(i4) :: h
integer(i4) :: hr

!-------------------------
! Compile met variables based on dayhour / nighthour

if (i == 1) then            ! Day time
  hr = dayhour
else if (i == 2) then       ! Night time
  hr = nighthour
end if

!-------------------------

Eg_b = 0.0
dEgdt_b = 0.0

! Calculate sensible heat flux between the air and soil surface (CLM eqn 5.60, CLM 4.0 Eqn 5.61)
! NOTE: Potential temperature of atmosphere substituted for simple surface air temperature

! if (i == 1) then          ! Day time
!
!   if (dayl > 0.) then
!
!     Datm = Patm / (Rspec * (tday + Tfreeze))
!
!     Hg_b    = -Datm * Cair * (tday + Tfreeze - Tsoil(snl+1)) / rah_b
!     dHgdt_b = Datm * Cair / rah_b
!
!     Swg_b = srad_dir * 1000. / (dayl * 3600.) / real(hr) / 10.
!     Lwg_b = lrad * 1000. / (dayl * 3600.) / real(hr) /10.
!
!     hs = Swg_b - Lwg_b - Hg_b
!     dhsdT = -dHgdt_b
!
!   else
!
!     Datm = Patm / (Rspec * (tday + Tfreeze))
!
!     Hg_b    = -Datm * Cair * (tday + Tfreeze - Tsoil(snl+1)) / rah_b
!     dHgdt_b = Datm * Cair / rah_b
!
!     Swg_b = 0.
!     Lwg_b = 0.
!
!     hs = 0.
!     dhsdT = -dHgdt_b
!
!   end if
!
! else if (i == 2) then     ! Night time
!
!   Datm = Patm / (Rspec * (tnight + Tfreeze))
!
!   Hg_b    = -Datm * Cair * (tnight + Tfreeze - Tsoil(snl+1)) / rah_b
!   dHgdt_b = Datm * Cair / rah_b
!
!   Swg_b = 0.0
!   Lwg_b = 0.0
!
!   hs = Swg_b - Lwg_b - Hg_b
!   hs = 0.
!   dhsdT = -dHgdt_b
!
! end if

!-------------------------
! Time-step dt = 3600 for now as subroutine written for constant 1 hour time-step
dt = 3600.

! TESTING - only one iteration per day night
dt = dt / 1.
hr = hr * 1
hs = hs / 1.

if (i == 1 .and. abs(tday-Tsoil(snl+1)+273.15) > 10.) then

  surf_atm_diff = abs(tday-Tsoil(snl+1)+273.15)

  dt = dt * 1. / surf_atm_diff
  hr = floor(hr / 1  * surf_atm_diff)
  hs = hs * 1. / surf_atm_diff

else if (i == 2 .and. abs(tnight-Tsoil(snl+1)+273.15) > 10.) then

  surf_atm_diff = abs(tnight-Tsoil(snl+1)+273.15)

  dt = dt * 1. / surf_atm_diff
  hr = floor(hr / 1  * surf_atm_diff)
  hs = hs * 1. / surf_atm_diff

end if

! if (i == 1) then
!
!   surf_atm_diff = abs(tday-Tsoil(snl+1)+273.15)
!
!   surf_atm_diff = max(min(10.,surf_atm_diff),1.)
!
!   dt = dt / surf_atm_diff
!   hr = ceiling(hr * surf_atm_diff)
!   hs = hs / surf_atm_diff
!
! else if (i == 2) then
!
!   surf_atm_diff = abs(tnight-Tsoil(snl+1)+273.15)
!
!   surf_atm_diff = max(min(10.,surf_atm_diff),1.)
!
!   dt = dt / surf_atm_diff
!   hr = ceiling(hr * surf_atm_diff)
!   hs = hs / surf_atm_diff
!
! end if

! print *, surf_atm_diff

!-------------------------
! Run water model in hourly time-step based on length of dayhour / nighthour
hourloop : do h = 1, hr

  ! Set-up tridiagonal system for current time-step soil temperature

  ! Top snow/soil layer is the layer-averaged temp, so to be more accurate
  ! we adjust the heat capacity of the top layer. This will give a more realistic top soil layer
  ! temperature. See P.107 of CLM 4.0. Eqn 6.29
  ! Start with top soil layer
  l = snl+1

  delZtop = 0.5 * (zpos(l) - zipos(l) + Tstune * (zpos(l+1) - zipos(l)))

  fact(l) = dt / Cl(l) * dz(l) / delZtop
  Fhti(l) = Kl(l) * (Tsoil(l+1) - Tsoil(l)) / (zpos(l+1) - zpos(l))

  ! Middle layers

  do l = snl+2, nl-1
    fact(l) = dt / Cl(l)
    Fhti(l) = -Kl(l) * (Tsoil(l) - Tsoil(l+1)) / (zpos(l+1) - zpos(l))
  end do

  ! Bottom soil layer

  fact(nl) = dt / Cl(nl)
  Fhti(nl) = 0.

  ! Set up the coefficients for the tridiagonal matrix solver
  ! top snow/soil layer
  l = snl+1

    dzp = zpos(l+1) - zpos(l)

    avect(l) = 0._dp
    bvect(l) = 1._dp + (1._dp - CNfac) * fact(l) * (Kl(l) / dzp) - fact(l) * dhsdT
    cvect(l) =       - (1._dp - CNfac) * fact(l) *  Kl(l) / dzp
    rvect(l) = Tsoil(l) + fact(l) * (hs - dhsdT * Tsoil(l) + CNfac * Fhti(l))

  ! Middle soil layers
  do l = snl+2, nl-1

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

    avect(l) =       - (1._dp - CNfac) * fact(l) *  Kl(l-1) / dzm
    bvect(l) = 1._dp + (1._dp - CNfac) * fact(l) *  Kl(l-1) / dzm
    cvect(l) = 0._dp
    rvect(l) = Tsoil(l) - CNfac * fact(l) * Fhti(l-1)

  ! Call tridiagonal solver
  a = snl+1
  b = nl

  call tridiag(avect(a:b),bvect(a:b),cvect(a:b),rvect(a:b),Tnsoi(a:b))

  ! print *, i, hs, hr, tday,tnight, Tsoil-Tfreeze, snl,zsno

  call soilphasechg(snl,Wsno,zsno,Wliq,Wice,Tsoil,Tnsoi,Tsat,Bexp,Psat,zpos,dz,dzmm,Kl,fact,Fhti,ithaw,&
                    Tice,Tpor,Tliq,fice0,fsnow,qsnomelt,hs,dhsdT,dt)

  Tsoiln(snl+1:nl) = Tsoil(snl+1:nl)
  Tsoil(snl+1:nl)  = Tnsoi(snl+1:nl)

end do hourloop

! print *, i, hs,hr, dayl, snl, tday,tnight, Tsoil(ns:nl)-Tfreeze, zsno, qsnomelt, dayl, srad_dir


end subroutine soiltemperature


end module soilphysicsmod
