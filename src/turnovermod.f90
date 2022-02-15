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

use pftparmod,    only : npft,pftpar
use statevarsmod, only : ndyear,dayvars,soilvars,vegvars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: small    = 1.

!-------------------------
! Pointer variables
logical,  pointer, dimension(:) :: present           ! PFT present
real(sp), pointer, dimension(:) :: nind              ! PFT population
real(sp), pointer, dimension(:) :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: turnover_ind      ! Total turnover of individual (gC m-2)

!-------------------------
! Local variables
real(sp) :: lm_turn         ! Leaf turnover of individual (gC m-2)
real(sp) :: rm_turn         ! Root turnover of individual (gC m-2)
real(sp) :: sm_turn         ! Sapwood turnover of individual (gC m-2)
real(sp) :: l_torate        ! Leaf turnover rate (yr-1)
real(sp) :: s_torate        ! Sapwood turnover rate (yr-1)
real(sp) :: r_torate        ! Root turnover rate (yr-1)
integer  :: pft

!-------------------------

present => vegvars(grid)%present
nind    => vegvars(grid)%nind
lm_ind  => vegvars(grid)%lm_ind
rm_ind  => vegvars(grid)%rm_ind
sm_ind  => vegvars(grid)%sm_ind
hm_ind  => vegvars(grid)%hm_ind

litter_ag_fast => vegvars(grid)%litter_ag_fast
litter_ag_slow => vegvars(grid)%litter_ag_slow
litter_bg      => vegvars(grid)%litter_bg
turnover_ind   => vegvars(grid)%turnover_ind

!-------------------------

do pft = 1, npft

  if (present(pft)) then

    l_torate = 1. / pftpar(8,pft)
    s_torate = 1. / pftpar(9,pft)
    r_torate = 1. / pftpar(10,pft)

    ! Calculate the biomass turnover in this year
    lm_turn = lm_ind(pft) * l_torate
    sm_turn = sm_ind(pft) * s_torate
    rm_turn = rm_ind(pft) * r_torate

    ! Update the carbon pools
    lm_ind(pft) = lm_ind(pft) - lm_turn
    sm_ind(pft) = sm_ind(pft) - sm_turn
    rm_ind(pft) = rm_ind(pft) - rm_turn

    ! Convert sapwood to heartwood
    hm_ind(pft) = hm_ind(pft) + sm_turn

    !Transfer to litter pools
    litter_ag_fast(pft) = litter_ag_fast(pft) + lm_turn * nind(pft)

    litter_bg(pft) = litter_bg(pft) + rm_turn * nind(pft)

    !Record total turnover of individual
    turnover_ind(pft) = lm_turn + sm_turn + rm_turn

  end if

end do

end subroutine turnover





end module turnovermod
