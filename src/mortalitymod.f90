module mortalitymod

! Module to calculate vegation mortality
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine mortality(mtemp,dtemp,abm_inc,sla,present,lm_ind,rm_ind,sm_ind,hm_ind,nind,litter_ag_fast,litter_ag_slow,&
                     litter_bg,turnover_ind)

use parametersmod, only : i4,sp
use pftparmod,     only : npft,pftpar,tree,boreal

real(sp), dimension(:), intent(in)    :: mtemp
real(sp), dimension(:), intent(in)    :: dtemp
real(sp), dimension(:), intent(in)    :: abm_inc
real(sp), dimension(:), intent(in)    :: sla               ! Specific leaf area (m2 gC-1)
logical,  dimension(:), intent(inout) :: present           ! PFT present
real(sp), dimension(:), intent(inout) :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: nind              ! PFT population
real(sp), dimension(:), intent(inout) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), dimension(:), intent(inout) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), dimension(:), intent(inout) :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), dimension(:), intent(inout) :: turnover_ind      ! Total turnover of individual (gC m-2)

!-------------------------
! Parameters
real(sp), parameter :: mort_max   = 0.01  ! Asymptotic maximum mortality rate (yr-1)
real(sp), parameter :: k_mort     = 0.3   ! Coefficient of growth efficiency in mortality equation
real(sp), parameter :: ramp_gddtw = 300.  ! Ramp for heat damage function

real(sp) :: gddtw
real(sp) :: mtemp_max
real(sp) :: twmax
real(sp) :: greffic
real(sp) :: bm_delta    ! Net individual living biomass increment (incorporating loss through leaf, root and sapwood turnover) (gC)
real(sp) :: mort        ! Tree mortality rate
real(sp) :: nind_kill   ! Reduction in individual density due to mortality (indiv/m2)
real(sp) :: litter_inc
real(sp) :: heatstress  ! Reduction in individual density (& establishment) due to heat induced mortality  (indiv/m2)

integer(i4) :: pft

!-------------------------

mtemp_max = maxval(mtemp)

!-------------------------

do pft = 1, npft

  twmax = pftpar%twmax(pft)

  if (present(pft) .and. tree(pft) .and. nind(pft) > 0.) then

    !Calculate net individual living biomass increment

    bm_delta = max(0., abm_inc(pft) / nind(pft) - turnover_ind(pft))

    !Calculate growth efficiency (net biomass increment per unit leaf area)

    greffic = bm_delta / lm_ind(pft) / sla(pft)

    !Mortality rate inversely related to growth efficiency (Prentice et al 1993)

    mort = mort_max / (1. + k_mort * greffic)

    !heat damage mortality in boreal trees

    if (mtemp_max > twmax) then  ! heat damage

      !calculate growing degree days above twmax

      gddtw = sum(dtemp - twmax, mask = dtemp > twmax)

      heatstress = min(1., gddtw / ramp_gddtw)

    else

      heatstress = 0.

    end if

    !------
    ! Reduce individual density (and thereby gridcell-level biomass) by mortality rate
    mort      = min(1., mort + heatstress)
    nind_kill = nind(pft) * mort
    nind(pft) = max(0., nind(pft) - nind_kill)

    !------
    ! Transfer lost biomass to litter
    litter_ag_fast(pft) = litter_ag_fast(pft) + nind_kill * lm_ind(pft)
    litter_ag_slow(pft) = litter_ag_slow(pft) + nind_kill *(sm_ind(pft) + hm_ind(pft))
    litter_bg(pft)      = litter_bg(pft)      + nind_kill * rm_ind(pft)

  end if

  !------

  if (nind(pft) == 0.) then

    present(pft)   = .false.
    lm_ind(pft)  = 0.
    sm_ind(pft)  = 0.
    hm_ind(pft)  = 0.
    rm_ind(pft)  = 0.

  end if

  !------

  ! write(0,*) 'PFT: ', pft, mort, heatstress, 'NIND_mort: ', nind_kill, nind(pft)

end do


end subroutine mortality


end module mortalitymod
