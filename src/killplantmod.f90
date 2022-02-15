module killplantmod

! Module to kill plant if productivity is too low
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine killplant(year,grid,day)

use pftparmod,    only : npft,pftpar,tree
use statevarsmod, only : ndyear,gppvars,vegvars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: bminc_ind_min = 1.e-6   ! Minimum annual productivity per individual has to be more than this value (gC)

!-------------------------
! Pointer variables
logical,  pointer, dimension(:) :: present           ! PFT present
real(sp), pointer, dimension(:) :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: nind              ! PFT population
real(sp), pointer, dimension(:) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: bm_inc

real(sp) :: bm_inc_ind
integer  :: pft

!-------------------------

present   => vegvars(grid)%present
lm_ind    => vegvars(grid)%lm_ind
rm_ind    => vegvars(grid)%rm_ind
sm_ind    => vegvars(grid)%sm_ind
hm_ind    => vegvars(grid)%hm_ind
nind      => vegvars(grid)%nind

litter_ag_fast => vegvars(grid)%litter_ag_fast
litter_ag_slow => vegvars(grid)%litter_ag_slow
litter_bg      => vegvars(grid)%litter_bg
bm_inc         => vegvars(grid)%bm_inc

!-------------------------


do pft = 1, npft

  if (present(pft)) then

    bm_inc_ind = bm_inc(pft) / nind(pft)

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
