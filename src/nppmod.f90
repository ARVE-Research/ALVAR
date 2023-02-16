module nppmod

! Module to calculate NPP coded for ALVAR
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)
! Implementation of carbon labile pool routine from Zaehle & Friend (2010)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine calcnpp(day,i,dayl,tday,tnight,tsoil,tree,present,gpp,dphen,leafon,leafondays,leafoffdays,lm_ind,rm_ind,sm_ind,hm_ind, &
                  nind,cellarea,areafrac,clabile,creserve,aresp,npp,npp_tot,bm_inc,abm_inc)

use parametersmod, only : i4,sp
use pftparmod,     only : npft,pftpar,raingreen,summergreen

implicit none

integer(i4),            intent(in)    :: day
integer(i4),            intent(in)    :: i
real(sp),               intent(in)    :: dayl         ! Daylength (h)
real(sp),               intent(in)    :: tday         ! Mean daytime temperature (degC)
real(sp),               intent(in)    :: tnight       ! Mean nighttime temperature (degC)
real(sp),               intent(in)    :: tsoil        ! Soil temperature (top layer) (K)
logical,  dimension(:), intent(in)    :: tree
logical,  dimension(:), intent(in)    :: present      ! PFT present
real(sp), dimension(:), intent(in)    :: gpp          ! Gross primary productivity under actual condition (g C m-2 d-1)
real(sp), dimension(:), intent(in)    :: dphen        ! Phenology status of summergreen (proportion of leaf-on) (fraction)
logical,  dimension(:), intent(in)    :: leafon       ! Leaf phenology on/off for photosynthesis
real(sp), dimension(:), intent(in)    :: leafondays   ! Number of days since leaf phenology is on
real(sp), dimension(:), intent(in)    :: leafoffdays    ! Number of days since leaf phenology is off
real(sp), dimension(:), intent(in)    :: lm_ind       ! Leaf carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(in)    :: rm_ind       ! Root carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(in)    :: sm_ind       ! Sapwood carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(in)    :: hm_ind       ! Sapwood carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(in)    :: nind         ! PFT population
real(sp),               intent(in)    :: cellarea     ! Area of gridcell (m2)
real(sp),               intent(in)    :: areafrac     ! Ground area fraction in gridcell (fraction)
real(sp), dimension(:), intent(inout) :: clabile      ! Carbon storage for growth respiration (gC m-2)
real(sp), dimension(:), intent(inout) :: creserve
real(sp), dimension(:), intent(inout) :: aresp        ! Autotrophic maintenence respiration (g C m-2 d-1)
real(sp), dimension(:), intent(inout) :: npp          ! Net primary productivity (g C m-2 d-1)
real(sp), dimension(:), intent(inout) :: npp_tot      ! Total net primary productivity (g C d-1)
real(sp), dimension(:), intent(inout) :: bm_inc
real(sp), dimension(:), intent(inout) :: abm_inc

!-------------------------
! Parameters
real(sp), parameter :: k       = 0.0548
real(sp), parameter :: tc      = 1. / 56.02   ! Constant in Arrhenius equation (Sitch et al. 2003; Eq. 23)
real(sp), parameter :: min_npp = 0.02
real(sp), parameter :: nlim    = 1.

!-------------------------
! Local variables
real(sp) :: temp
real(sp) :: tsoil_C               ! Soil temperature (degC)
real(sp) :: gtemp_air             ! Arrhenius equation varaible for air
real(sp) :: gtemp_soil            ! Arrhenius equation varaible for soil
real(sp) :: lresp_ind             ! Leaf respiration for individual (gC m-2 d-1)
real(sp) :: rresp_ind             ! Root respiration for individual (gC m-2 d-1)
real(sp) :: sresp_ind             ! Sapwood respiration for individual (gC m-2 d-1)
real(sp) :: aresp_ind             ! Total respiration for individual (gC m-2 d-1)

real(sp) :: diurfrac
real(sp) :: respcoeff
real(sp) :: l_c2n
real(sp) :: r_c2n
real(sp) :: s_c2n
real(sp) :: crve_time

real(sp) :: dphen_g
real(sp) :: biomass_tot

integer(i4) :: pft

!-------------------------

if (areafrac <= 0.) then
  npp = -9999.
  npp_tot = 0.
  return
end if

!-------------------------

if (i == 1) then
  diurfrac = dayl / 24.
  temp = tday
else
  diurfrac = 1. - (dayl / 24.)
  temp = tnight
end if

!-------------------------

tsoil_C = tsoil - Tfreeze

! Calculate variables for Arrhenius equation (Sitch et al., 2003; Eq. 23)
if (temp > -40.) then
  gtemp_air = exp(308.56 * (tc - 1. / (temp  + 46.02)))
else
  gtemp_air = 0.
end if

if (tsoil_C > -40.) then
  gtemp_soil = exp(308.56 * (tc - 1. / (tsoil_C  + 46.02)))
else
  gtemp_soil = 0.
end if

!-------------------------

do pft = 1, npft

  if (present(pft)) then

    !------

    respcoeff = pftpar%respcoeff(pft)
    l_c2n     = pftpar%l_c2n(pft)
    r_c2n     = pftpar%r_c2n(pft)

    if (tree(pft)) then
      s_c2n = 1. / pftpar%s_c2n(pft)   ! To avoid divide by zero in grass PFTs
    else
      s_c2n = 0.
    end if

    !------

    ! Calculate respiratoin rate for each carbon compartment
    ! Multiply daylength fraction (diurfrac) to convert day-1 rate into day/night
    lresp_ind = (lm_ind(pft) / l_c2n) * respcoeff * k * gtemp_air  * dphen(pft) * diurfrac
    rresp_ind = (rm_ind(pft) / r_c2n) * respcoeff * k * gtemp_soil * dphen(pft) * diurfrac
    sresp_ind = (sm_ind(pft) * s_c2n) * respcoeff * k * gtemp_air  * diurfrac

    ! Sum all maintenence respiration rate of individual
    aresp_ind = lresp_ind + rresp_ind + sresp_ind

    ! Multiply individual respiration with population to find total repiration of the day
    aresp(pft) = aresp_ind * nind(pft)

    ! if (isnan(aresp(pft))) write(0,*)'NaN respiration at',gridlon(grid),gridlat(grid),dayl,gtemp_air,diurfrac,aresp(pft),temp

    !-------------------------
    ! Subtract respiration from GPP to get NPP (gC m-2 d-1)
    ! NOTE: Included new carbon reserve routine for summergreen and herbaceous PFTs

    if (i == 1) then  ! Day

      npp(pft) = gpp(pft) - aresp(pft)

    else  ! Night

      npp(pft) = npp(pft) - aresp(pft)

      !-------------------------
      ! Make up negative NPP with carbon storage for growht respirtaion
      if (summergreen(pft) .or. .not.tree(pft)) then

        crve_time = pftpar%crve_time(pft)     ! Get the duration of onset allowed for access creserve (based on ORCHIDEE) (Krinner et al., 2005)

        if (leafon(pft)) then

          !-------------------------
          ! First day of leaf onset
          if (leafondays(pft) == 1) then

            ! Save proportion of reserve into labile for translocation in the coming days
            clabile(pft)  = creserve(pft) / crve_time

            creserve(pft) = creserve(pft) - clabile(pft)
            npp(pft)      = npp(pft) + clabile(pft)

          !-------------------------
          ! Day 1 to 30 of leaf onset, allow access to carbon reserve from previous year
          else if (leafondays(pft) > 1 .and. leafondays(pft) <= crve_time) then

            creserve(pft) = creserve(pft) - clabile(pft)
            npp(pft)      = npp(pft) + clabile(pft)

            if (leafondays(pft) == crve_time) clabile(pft) = 0.

          !-------------------------
          ! Day 31 and beyond, accumulate NPP back into carbon reserve
          else if (leafondays(pft) > crve_time) then

            if (npp(pft) > 0.) then

              creserve(pft) = creserve(pft) + 0.25 * npp(pft)
              creserve(pft) = max(creserve(pft), 0.)

              npp(pft)     = 0.75 * npp(pft)

            end if ! NPP > 0

          end if ! Leadondays

        end if ! Leafon

        !-------------------------
        ! End of growing season, approximate reserve lost over the winter / unfavorbale conditions
        if (summergreen(pft) .and. leafoffdays(pft) == 1) clabile(pft) = clabile(pft) * 0.7
        if (.not.tree(pft)   .and. leafoffdays(pft) >= 1) clabile(pft) = clabile(pft) * 0.9981

        ! if (leafoffdays(pft) >= 1) clabile(pft) = clabile(pft) * 0.9981

        bm_inc(pft) = npp(pft)

      else    ! Other evergreen tree PFTs

        bm_inc(pft) = npp(pft)

      end if  ! Summergreen / herbaceous

      !------

      if (day == 1) then
        abm_inc(pft) = 0.
        abm_inc(pft) = abm_inc(pft) + bm_inc(pft)
      else
        abm_inc(pft) = abm_inc(pft) + bm_inc(pft)
      end if

      !-------------------------
      ! Calculate total daily NPP of entire gridcell (gC d-1)
      npp_tot(pft) = npp(pft) * cellarea * areafrac

      biomass_tot = lm_ind(pft) + rm_ind(pft) + sm_ind(pft)

      !-------------------------

      if (areafrac <= 0.) then
        npp(pft) = -9999
        npp_tot(pft) = 0.
      end if

    end if ! Day/night IF loop

  end if ! Present IF loop

end do ! PFT loop


end subroutine calcnpp

!---------------------------------------------------------------------

end module nppmod
