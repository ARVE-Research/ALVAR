module soilphasechgmod

! Module adapted from ARVE-DGVM to calculate soil \phase change (Leo Lai, Aug 2021)

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine soilphasechg(snl,Wsno,zsno,Wliq,Wice,Tsoil,Tnsoi,Tsat,Bexp,Psat,zpos,dz,dzmm,Kl,fact,Fhti,ithaw,&
                        Tice,Tpor,Tliq,fice0,fsnow,qsnomelt,hs,dhsdT,dt)
!calculates the energy budget of phase change for the soil water constituents

! use soilstate_vars, only : hs,dhsdt

use parametersmod, only : i4,sp,dp,Tfreeze
use statevarsmod, only : nl,ns

implicit none

integer(i4),                  intent(in)    :: snl     !index value of top snow layer (negative is more layers)
real(sp),                     intent(inout) :: Wsno   !snow water equivalent of the snowpack (mm)
real(sp),                     intent(inout) :: zsno   !total thickness of the snowpack (m)
real(sp),   dimension(ns:nl), intent(inout) :: Wliq    !soil liquid water content at layer midpoint (mm)
real(sp),   dimension(ns:nl), intent(inout) :: Wice    !soil ice content at layer midpoint (mm)
real(sp),   dimension(ns:nl), intent(inout) :: Tsoil   !soil temperature (K)
real(dp),   dimension(ns:nl), intent(inout) :: Tnsoi   !updated soil temperature at timestep n+1
real(sp),   dimension(ns:nl), intent(inout) :: Tsat    !soil water volumetric water content at saturation (fraction)
real(sp),   dimension(ns:nl), intent(inout) :: Bexp    !soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(sp),   dimension(ns:nl), intent(inout) :: Psat    !soil water matric potential at saturation (mm)
real(sp),   dimension(ns:nl), intent(in)    :: zpos    !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(sp),   dimension(ns:nl), intent(in)    :: dz      !thickness of the soil layers (m)
real(sp),   dimension(ns:nl), intent(in)    :: dzmm    !thickness of the soil layers (mm)
real(sp),   dimension(ns:nl), intent(inout) :: Kl         ! Soil thermal conductivity across layer boundary (W m-1 K-1)
real(sp),   dimension(ns:nl), intent(inout) :: fact   !factor used in computing tridiagonal coefficients
real(sp),   dimension(ns:nl), intent(inout) :: Fhti     !heat flux across soil layer boundary (W m-2)
integer(i4), dimension(ns:nl), intent(inout) :: ithaw
real(sp),   dimension(ns:nl), intent(inout) :: Tice           !soil ice content (fraction)
real(sp),   dimension(ns:nl), intent(inout) :: Tpor    !soil volumetric porosity (fraction)
real(sp),   dimension(ns:nl), intent(inout) :: Tliq    !soil liquid water content (fraction)
real(sp),   dimension(ns:nl), intent(inout) :: fice0   !layer ice fraction, previous timestep
real(sp),                     intent(inout) :: fsnow   !fraction of the gridcell covered by snow (fraction)
real(sp),                     intent(inout) :: qsnomelt ! Snow melt (kg m-2 s-1)
real(sp),                     intent(inout) :: hs       ! Net energy flux into the surface (W m-2)
real(sp),                     intent(inout) :: dhsdT
real(sp),                     intent(inout) :: dt

!parameters
real(sp), parameter :: psno_new        = 100.        ! density of new snow (kg m-3)
integer, parameter :: m = 1                             ! exponent for fsnow calc. Suggested to be 1 for global applications
real(sp), parameter :: Lf    = 3.337e5          !Water latent heat of fusion (J kg-1)
real(sp), parameter :: grav  = 9.80616       !Gravitational acceleration  (m s-2)
real(sp), parameter :: CNfac  = 0.5          !Crank Nicholson factor between 0 and 1
real(sp), parameter :: oneminusCNfac = 1. - CNfac
real(sp), parameter :: pliq  = 1000.         !density of water (kg m-3)
real(sp), parameter :: pice  =  917.         !density of ice (kg m-3)
real(sp), parameter :: z0mg  = 0.01          !momentum roughness length for bare soil (m, CLM eqn 3.49)



!local variables
real(sp), dimension(ns:nl) :: Hi     !excess heat before phase change calculation (W m-2)
real(sp), dimension(ns:nl) :: Hin    !excess heat after phase change calculation (W m-2)
real(sp), dimension(ns:nl) :: Wliqn  !quantity of liquid water in the soil at time t+1, after phase change calc. (mm)
real(sp), dimension(ns:nl) :: Wicen  !quantity of ice in the soil at time t+1, after phase change calc. (mm)
real(sp), dimension(ns:nl) :: Wtot   !total mass of liquid and ice in the soil (mm)
real(sp), dimension(ns:nl) :: thaw   !quantity of water thawed, negative means frozen (mm)
real(sp), dimension(ns:nl) :: Fnhti  !heat flux across soil layer boundary at time t+1 (W m-2)
real(sp), dimension(ns:nl) :: Bi     !temporary variable for phase change calculations
real(sp) :: Ephase    !total latent heat of phase change W m-2
real(sp) :: dTemp     !change in temperature of the soil layer between timesteps
real(sp) :: frac
real(sp) :: tmp1
real(sp) :: Tliqmax
real(sp), dimension(nl) :: Wliqmax   !maximum liquid water fraction when the soil temp is below freezing   (mm)
integer(i4) :: l            !counters

!----------------------------
!calculate soil energy excess or deficit vs. water phase change

!fn1 in CLM (6.24)
do l = snl+1,nl-1
    Fnhti(l) =  Kl(l) * (Tnsoi(l+1) - Tnsoi(l)) / (zpos(l+1) - zpos(l))   !CLM 6.24
end do

Fnhti(nl) = 0.  !no heat flux in bottom layer

!Brr in CLM  (middle two terms of CLM 6.42)
   !top snow/soil layer
    Bi(snl+1) = CNfac * Fhti(snl+1) + oneminusCNfac * Fnhti(snl+1)

    !other layers
do l = snl+2,nl
    Bi(l) = CNfac * (Fhti(l) - Fhti(l-1)) + oneminusCNfac * (Fnhti(l) - Fnhti(l-1))
end do

        !Frozen soil method by Niu & Yang, 2006, J Hydrometeorology (7) 937-952.
        !coded by Joe Melton Dec 14 2007

do l = 1,nl !soil layers only
 if (Tnsoi(l) < Tfreeze) then

  !(Niu & Yang (eqn 3)

    Tliqmax = Tsat(l) * (1.e3 * Lf * (Tnsoi(l) - Tfreeze) / (grav * Psat(l) * Tnsoi(l)))**(-1. / Bexp(l))

    if (Tliqmax > Tsat(l)) Tliqmax = Tsat(l)

    Wliqmax(l) = Tliqmax * dzmm(l)

 end if
end do

!--

Wliqn(:) = Wliq(:)
Wicen(:) = Wice(:)

Ephase   = 0.

!--

do l = snl+1,nl  !FLAG, not sure about whether this should be only gnl or not...

    ithaw(l) = 0
    Hi(l)    = 0.
    thaw(l)  = 0.
    Wtot(l)  = Wice(l) + Wliq(l)

           if (Wtot(l) > 0.) then
             fice0(l) = Wice(l) / Wtot(l)
           else  !no water/ice so skip
             fice0(l) = 0.
             cycle
           end if

  !assess conditions
          if (l >= 1) then ! soil layers

           if (Tnsoi(l) > Tfreeze .and. Wice(l) > 0.) then       !thawing conditions
               ithaw(l) = 1
               Tnsoi(l) = Tfreeze
             else if (Tnsoi(l) < Tfreeze .and. Wliq(l) > Wliqmax(l)) then  !freezing conditions
               ithaw(l) = 2
               Tnsoi(l) = Tfreeze
             else                                                 !neither freeze nor thaw
               ithaw(l) = 0
           end if

           if (snl+1 == 1 .and. Wsno > 0. .and. l == 1) then     !conditions where there is snow but no explict snow layer
             if (Tsoil(l) > Tfreeze) then
               ithaw = 1  !thaw
               Tnsoi(l) = Tfreeze
             end if
           end if

         else !snow layers
           if (Tnsoi(l) > Tfreeze .and. Wice(l) > 0.) then       !thawing conditions
               ithaw(l) = 1
               Tnsoi(l) = Tfreeze
             else if (Tnsoi(l) < Tfreeze .and. Wliq(l) > 0.) then  !freezing conditions
               ithaw(l) = 2
               Tnsoi(l) = Tfreeze
             else                                                 !neither freeze nor thaw
               ithaw(l) = 0
           end if

           if (snl+1 == 1 .and. Wsno > 0. .and. l == 1) then     !conditions where there is snow but no explict snow layer
             if (Tsoil(l) > Tfreeze) then
               ithaw = 1
               Tnsoi(l) = Tfreeze
             end if
           end if


         end if

  !calculate amount of water frozen or thawed
           if (ithaw(l) /= 0) then   !if there is freeze or thaw

                 dTemp = Tnsoi(l) - Tsoil(l)

                 !calculate heat excess (deficit)   (CLM 6.42)
                 if (l > snl+1) then
                  !all layers except surface
                   Hi(l) = Bi(l) - dTemp / fact(l)
                 else
                   !surface layer
                   Hi(l) = hs + dhsdT * dTemp + Bi(l) - dTemp / fact(l)
                 end if

           end if

           !if there are thawing (freezing) conditions but the heat excess (deficit) is negative (positive),
           !then there is no thaw (freeze)

           !Notes: Hi(l) = CLM hm, thaw(l) = CLM xm, xmf = Ephase

           if ((ithaw(l) == 1 .and. Hi(l) < 0.) .or. (ithaw(l) == 2 .and. Hi(l) > 0.)) then
             Hi(l) = 0.
             ithaw(l) = 0
           end if

           if (ithaw(l) /= 0 .and. abs(Hi(l)) > 0.) then


                thaw(l) = Hi(l) * dt / Lf           !amount of water thawed (frozen) kg m-2

                !----
                !if snow exists but there is no explicit layer

                if (l == 1) then
                  if (snl+1 == 1 .and. Wsno > 0. .and. thaw(l) > 0.) then

                    tmp1 = Wsno
                    Wsno = max(0.,tmp1 - thaw(l))
                    frac = Wsno / tmp1
                    zsno = frac * zsno

                    ! fractional area of the gridcell covered by snow (Niu and Yang 2007)
                    if (Wsno > 0.) then
                      fsnow = tanh(zsno / (2.5_dp * z0mg * ((min(Wsno/zsno, 800.)) / psno_new)**m))
                    else
                      fsnow = 0.
                    end if

                    Hin(l) = Hi(l) - Lf * (tmp1 - Wsno) / dt

                    if (Hin(l) > 0.) then
                      thaw(l) = Hin(l) * dt / Lf
                      Hi(l)   = Hin(l)
                    else
                      thaw(l) = 0.
                      Hi(l)   = 0.
                    end if

                    qsnomelt = max(0., (tmp1 - Wsno))/ dt  !kg m-2 s-1
                    Ephase = Lf * qsnomelt
                  end if
                end if

    !calculate additional heat excess (deficit) after thawing (freezing)
             if (l > 0) then !soil layers only

                  Hin(l) = 0.

                      if (thaw(l) > 0.) then  !thawing (CLM 6.43)

                             Wicen(l) = max(0.,(Wice(l) - thaw(l)))

                     Hin(l) = Hi(l) - Lf * (Wice(l) - Wicen(l)) / dt


                   else if (thaw(l) < 0.) then !freezing  (CLM 6.43)

                            Wicen(l) = min(Wtot(l),Wice(l) - thaw(l))

                     Hin(l) = Hi(l) - Lf * (Wice(l) - Wicen(l)) / dt

                   end if


            else  !snow!

                  Hin(l) = 0.

                  if (thaw(l) > 0.) then  !thawing

                    Wicen(l) = max(0., Wice(l) - thaw(l))
                    Hin(l) = Hi(l) - Lf * (Wice(l) - Wicen(l)) / dt

                  else if (thaw(l) < 0.) then  !freezing

                    Wicen(l) = min(Wtot(l), Wice(l) - thaw(l))
                    Hin(l) = Hi(l) - Lf * (Wice(l) - Wicen(l)) / dt

                  end if

                end if

      Wliqn(l) = max(0., Wtot(l) - Wicen(l))

        !----------------------------------------------------------------

                if (abs(Hin(l)) > 0.) then

                  !if there is excess (deficit) heat, update the soil temperature

                  if (l > snl+1) then
                    Tnsoi(l) = Tnsoi(l) + fact(l) * Hin(l)

                  else
                    Tnsoi(l) = Tnsoi(l) + fact(l) * Hin(l) / (1. - fact(l) * dhsdT)

                  end if

                  !if there are both ice and water in the soil the temperature must be 0 deg C
                   !=---NOTE: This constraint is now removed with the Niu & Yang correction implemented above! -JM Dec 14 2007
                  !if (Wliqn(l) * Wicen(l) > 0.) Tnsoi(l) = Tfreeze

                end if

                Ephase = Ephase + Lf * (Wice(l) - Wicen(l)) / dt

                if (l < 1 .and. ithaw(l) == 1) then
                  qsnomelt = max(0., (Wice(l) - Wicen(l))) / dt  !FLAG does this get used anywhere?
                end if

           end if


end do

!-----
Wliq(:) = Wliqn(:)
Wice(:) = Wicen(:)

  !update the soil volumetric ice content, porosity and water content.
  Tice(1:nl) = min(Tsat(1:nl), Wice(1:nl) / (dz(1:nl) * pice))
  Tpor(1:nl) = Tsat(1:nl) - Tice(1:nl)
  Tpor(1:nl) = max(0.,Tpor(1:nl))
  Tliq(1:nl) = min(Tpor(1:nl), Wliq(1:nl) / (dz(1:nl) * pliq))

return

end subroutine soilphasechg

!---------------------------------------------------------------------

end module soilphasechgmod
