module mortalitymod

! Module to calculate vegation mortality
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine mortality(year,grid,day)

use pftparmod,    only : npft,pftpar,tree,boreal
use statevarsmod, only : ndyear,genvars,dayvars,soilvars,gppvars,vegvars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: mort_max   = 0.01  ! Asymptotic maximum mortality rate (yr-1)
real(sp), parameter :: k_mort     = 0.3   ! Coefficient of growth efficiency in mortality equation
real(sp), parameter :: ramp_gddtw = 300.  ! Ramp for heat damage function

!-------------------------
! Pointer variables
logical,  pointer, dimension(:) :: present           ! PFT present
real(sp), pointer, dimension(:) :: sla               ! Specific leaf area (m2 gC-1)
real(sp), pointer, dimension(:) :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: nind              ! PFT population
real(sp), pointer, dimension(:) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: turnover_ind      ! Total turnover of individual (gC m-2)
real(sp), pointer, dimension(:) :: bm_inc

real(sp),              dimension(12) :: mtemp
real(sp), allocatable, dimension(:)  :: dtemp

real(sp) :: gddtw
real(sp) :: mtemp_max
real(sp) :: twmax
real(sp) :: greffic
real(sp) :: bm_delta    ! Net individual living biomass increment (incorporating loss through leaf, root and sapwood turnover) (gC)
real(sp) :: mort        ! Tree mortality rate
real(sp) :: nind_kill   ! Reduction in individual density due to mortality (indiv/m2)
real(sp) :: litter_inc
real(sp) :: heatstress  ! Reduction in individual density (& establishment) due to heat induced mortality  (indiv/m2)

integer :: pft

!-------------------------

present   => vegvars(grid)%present
sla       => vegvars(grid)%sla
lm_ind    => vegvars(grid)%lm_ind
rm_ind    => vegvars(grid)%rm_ind
sm_ind    => vegvars(grid)%sm_ind
hm_ind    => vegvars(grid)%hm_ind
nind      => vegvars(grid)%nind

litter_ag_fast => vegvars(grid)%litter_ag_fast
litter_ag_slow => vegvars(grid)%litter_ag_slow
litter_bg      => vegvars(grid)%litter_bg
turnover_ind   => vegvars(grid)%turnover_ind
bm_inc         => vegvars(grid)%bm_inc

!-------------------------

allocate(dtemp(ndyear))

mtemp = genvars%tmp(5:16)
dtemp = dayvars(grid,1:ndyear)%tmean

mtemp_max = maxval(mtemp)

!-------------------------

do pft = 1, npft

  twmax = pftpar(30,pft)

  if (present(pft) .and. tree(pft) .and. nind(pft) > 0.) then

    !Calculate net individual living biomass increment

    bm_delta = max(0., bm_inc(pft) / nind(pft) - turnover_ind(pft))

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

  if (lprint .and. grid==gprint) write(0,*) 'PFT: ', pft, mort, heatstress, &
                          'NIND_mort: ', nind_kill, nind(pft), vegvars(grid)%fpc_grid(pft)

end do


end subroutine mortality


end module mortalitymod
