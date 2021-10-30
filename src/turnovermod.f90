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
real(sp), parameter :: small = 1.

!-------------------------
! Pointer variables
logical, pointer :: present
real(sp), pointer :: nind
real(sp), pointer :: hm_ind
real(sp), pointer :: lm_ind
real(sp), pointer :: sm_ind
real(sp), pointer :: rm_ind
real(sp), pointer :: litter_ag_fast
real(sp), pointer :: litter_ag_slow
real(sp), pointer :: litter_bg
real(sp), pointer :: turnover_ind

!-------------------------
! Local variables
! integer  :: pft
real(sp) :: l_torate
real(sp) :: s_torate
real(sp) :: r_torate
real(sp) :: lm_turn
real(sp) :: sm_turn
real(sp) :: rm_turn

!-------------------------

present => vegvars(grid,day)%present
nind => vegvars(grid,day)%nind
hm_ind => vegvars(grid,day)%hm_ind
lm_ind => vegvars(grid,day)%lm_ind
sm_ind => vegvars(grid,day)%sm_ind
rm_ind => vegvars(grid,day)%rm_ind

litter_ag_fast => vegvars(grid,day)%litter_ag_fast
litter_ag_slow => vegvars(grid,day)%litter_ag_slow
litter_bg => vegvars(grid,day)%litter_bg
turnover_ind => vegvars(grid,day)%turnover_ind

!-------------------------

if (present) then

  l_torate = 1.         ! LPJ turnover rate based on tissue logevity (Sitch et al., 2003 eq.26)
  r_torate = 1.
  s_torate = 0.05

  !Calculate the biomass turnover in this year

  lm_turn = lm_ind * l_torate
  sm_turn = sm_ind * s_torate
  rm_turn = rm_ind * r_torate

  !Update the pools

  lm_ind = lm_ind - lm_turn
  sm_ind = sm_ind - sm_turn
  rm_ind = rm_ind - rm_turn

  !Convert sapwood to heartwood

  hm_ind = hm_ind + sm_turn

  !Transfer to litter pools

  litter_ag_fast = litter_ag_fast + lm_turn * nind

  litter_bg = litter_bg + rm_turn * nind

  !Record total turnover

  turnover_ind = lm_turn + sm_turn + rm_turn

end if

end subroutine turnover





end module turnovermod
