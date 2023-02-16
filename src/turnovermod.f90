module turnovermod

! Module to calculate vegetation carbon turnover
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine turnover(present,nind,lm_ind,rm_ind,sm_ind,hm_ind,litter_ag_fast,litter_ag_slow,litter_bg,turnover_ind)

use parametersmod, only : i4,sp
use pftparmod,     only : npft,pftpar

implicit none

logical,  dimension(:), intent(in)    :: present           ! PFT present
real(sp), dimension(:), intent(in)    :: nind              ! PFT population
real(sp), dimension(:), intent(inout) :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), dimension(:), intent(inout) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), dimension(:), intent(inout) :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), dimension(:), intent(inout) :: turnover_ind      ! Total turnover of individual (gC m-2)

!-------------------------
! Parameters
real(sp), parameter :: small = 1.

!-------------------------
! Local variables
real(sp)    :: lm_turn         ! Leaf turnover of individual (gC m-2)
real(sp)    :: rm_turn         ! Root turnover of individual (gC m-2)
real(sp)    :: sm_turn         ! Sapwood turnover of individual (gC m-2)
real(sp)    :: l_torate        ! Leaf turnover rate (yr-1)
real(sp)    :: s_torate        ! Sapwood turnover rate (yr-1)
real(sp)    :: r_torate        ! Root turnover rate (yr-1)
integer(i4) :: pft

!-------------------------

do pft = 1, npft

  if (present(pft)) then

    l_torate = 1. / pftpar%l_torate(pft)
    s_torate = 1. / pftpar%s_torate(pft)
    r_torate = 1. / pftpar%r_torate(pft)

    ! Calculate the biomass turnover in this year
    lm_turn = lm_ind(pft) * (l_torate / 365.)
    sm_turn = sm_ind(pft) * (s_torate / 365.)
    rm_turn = rm_ind(pft) * (r_torate / 365.)

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
