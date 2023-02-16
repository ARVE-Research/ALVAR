module killplantmod

! Module to kill plant if productivity is too low
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine killplant(abm_inc,present,lm_ind,rm_ind,sm_ind,hm_ind,nind,litter_ag_fast,litter_ag_slow,litter_bg)

use parametersmod, only : i4,sp
use pftparmod,     only : npft,tree

!-------------------------
! Pointer variables
real(sp), dimension(:), intent(in)    :: abm_inc
logical,  dimension(:), intent(inout) :: present           ! PFT present
real(sp), dimension(:), intent(inout) :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), dimension(:), intent(inout) :: nind              ! PFT population
real(sp), dimension(:), intent(inout) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), dimension(:), intent(inout) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), dimension(:), intent(inout) :: litter_bg         ! Below ground litter pool (gC m-2)

!-------------------------
! Parameters
real(sp), parameter :: bminc_ind_min = 1.e-6   ! Minimum annual productivity per individual has to be more than this value (gC)

real(sp)    :: bm_inc_ind
integer(i4) :: pft

!-------------------------

do pft = 1, npft

  if (present(pft)) then

    ! bm_inc_ind = bm_inc(pft) / nind(pft)

    bm_inc_ind = abm_inc(pft) / nind(pft) ! Get annual biomass increment (abm_inc) as in LPJ

    if (bm_inc_ind < bminc_ind_min) then

      ! Not enough C increment this year, kill PFT and transfer carbon to litter

      ! All PFTs leaf and root mass
      litter_ag_fast(pft) = litter_ag_fast(pft) + lm_ind(pft) * nind(pft)
      litter_bg(pft)      = litter_bg(pft)      + rm_ind(pft) * nind(pft)

      lm_ind(pft) = 0.
      rm_ind(pft) = 0.

      if (tree(pft)) then  ! Remove the stem mass for woody PFTs

        litter_ag_slow(pft) = litter_ag_slow(pft) + (sm_ind(pft) + hm_ind(pft)) * nind(pft)

        sm_ind(pft) = 0.
        hm_ind(pft) = 0.

      end if

      ! Reset nind and present
      nind(pft)    = 0.
      present(pft) =.false.
      lm_ind(pft) = 0.
      sm_ind(pft) = 0.
      hm_ind(pft) = 0.
      rm_ind(pft) = 0.

    end if

  end if

end do



end subroutine killplant


end module killplantmod
