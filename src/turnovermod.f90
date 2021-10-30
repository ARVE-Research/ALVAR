module turnovermod

! Module to calculate vegetation carbon turnover
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine turnover(year,grid,day)

use statevarsmod, only : ndyear,dayvars,soilvars,vegvars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: small    = 1.
real(sp), parameter :: l_torate = 0.5      ! Leaf turnover rate (yr-1) NOTE: LPJ turnover rate based on tissue logevity (Sitch et al., 2003; Eq. 26)
real(sp), parameter :: r_torate = 0.5      ! Root turnover rate (yr-1)
real(sp), parameter :: s_torate = 0.05    ! Sapwood turnover rate (yr-1)

!-------------------------
! Pointer variables
logical,  pointer :: present           ! PFT present
real(sp), pointer :: nind              ! PFT population
real(sp), pointer :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), pointer :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), pointer :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), pointer :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), pointer :: turnover_ind      ! Total turnover of individual (gC m-2)

!-------------------------
! Local variables
real(sp) :: lm_turn         ! Leaf turnover of individual (gC m-2)
real(sp) :: rm_turn         ! Root turnover of individual (gC m-2)
real(sp) :: sm_turn         ! Sapwood turnover of individual (gC m-2)
! integer  :: pft
! real(sp) :: l_torate
! real(sp) :: s_torate
! real(sp) :: r_torate

!-------------------------

present => vegvars(grid,day)%present
nind    => vegvars(grid,day)%nind
lm_ind  => vegvars(grid,day)%lm_ind
rm_ind  => vegvars(grid,day)%rm_ind
sm_ind  => vegvars(grid,day)%sm_ind
hm_ind  => vegvars(grid,day)%hm_ind

litter_ag_fast => vegvars(grid,day)%litter_ag_fast
litter_ag_slow => vegvars(grid,day)%litter_ag_slow
litter_bg      => vegvars(grid,day)%litter_bg
turnover_ind   => vegvars(grid,day)%turnover_ind

!-------------------------

if (present) then

  ! Calculate the biomass turnover in this year
  lm_turn = lm_ind * l_torate
  sm_turn = sm_ind * s_torate
  rm_turn = rm_ind * r_torate

  ! Update the carbon pools
  lm_ind = lm_ind - lm_turn
  sm_ind = sm_ind - sm_turn
  rm_ind = rm_ind - rm_turn

  ! Convert sapwood to heartwood
  hm_ind = hm_ind + sm_turn

  !Transfer to litter pools
  litter_ag_fast = litter_ag_fast + lm_turn * nind

  litter_bg = litter_bg + rm_turn * nind

  !Record total turnover of individual
  turnover_ind = lm_turn + sm_turn + rm_turn

end if

end subroutine turnover





end module turnovermod
