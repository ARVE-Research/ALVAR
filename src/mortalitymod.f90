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

use statevarsmod, only : ndyear,dayvars,soilvars,vegvars,lprint,gprint

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
logical,  pointer :: present           ! PFT present
real(sp), pointer :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), pointer :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer :: nind              ! PFT population
real(sp), pointer :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), pointer :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), pointer :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), pointer :: turnover_ind      ! Total turnover of individual (gC m-2)

real(sp) :: bm_inc_ind

!-------------------------

lm_ind    => vegvars(grid,day)%lm_ind
rm_ind    => vegvars(grid,day)%rm_ind
sm_ind    => vegvars(grid,day)%sm_ind
hm_ind    => vegvars(grid,day)%hm_ind
nind      => vegvars(grid,day)%nind

litter_ag_fast => vegvars(grid,day)%litter_ag_fast
litter_ag_slow => vegvars(grid,day)%litter_ag_slow
litter_bg      => vegvars(grid,day)%litter_bg
turnover_ind   => vegvars(grid,day)%turnover_ind

!-------------------------


end subroutine mortality


end module mortalitymod
