module lightmod

! Module to calculate light competition aong vegetation and constrain grid FPC
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine light(year,grid,day)

use pftparmod,    only : npft,pftpar,tree
use statevarsmod, only : vegvars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: fpc_tree_max = 0.95    ! Maximum tree FPC (Sitch et al., 2003; P. 170)

!-------------------------
! Pointer variables
logical,  pointer, dimension(:) :: present           ! PFT present
real(sp), pointer, dimension(:) :: fpc_grid          ! Foilage projective cover over grid (fraction)
real(sp), pointer, dimension(:) :: fpc_ind           ! Foliage projective cover of individual (fraction)
real(sp), pointer, dimension(:) :: fpc_inc           ! Foliage projective cover increment (fraction)
real(sp), pointer, dimension(:) :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: nind              ! PFT population
real(sp), pointer, dimension(:) :: sla               ! Specific leaf area (m2 gC-1)
real(sp), pointer, dimension(:) :: crownarea         ! Tree crownarea (m2)
real(sp), pointer, dimension(:) :: lai_ind           ! Leaf area index of individual (m2 m-2)
real(sp), pointer, dimension(:) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: turnover_ind      ! Total turnover of individual (gC m-2)

!-------------------------
! Local variables
integer :: pft
integer :: ntree             ! Number of tree PFTs currently present
integer :: ngrass            ! Number of tree PFTs currently present

real(sp) :: fpc_inc_tree     ! Current years total FPC increment for tree PFTs
real(sp) :: fpc_tree_total   ! Total grid FPC for tree PFTs
real(sp) :: fpc_grass_max    ! Max allowed grass FPC given the current tree cover
real(sp) :: fpc_grass_total  ! Total grid FPC for grass PFTs
real(sp) :: grasscover       ! Grass PFT proportional cover ("crownarea" equivalent for grass)
real(sp) :: excess           ! Total tree FPC or grass cover to be reduced

real(sp) :: nind_kill        ! Reduction in individual density to reduce tree FPC to permitted maximum (indiv m-2)
real(sp) :: rm_kill          ! Reduction in grass PFT root mass to reduce grass cover to permitted maximum (gC)
real(sp) :: lm_kill          ! Reduction in grass PFT leaf mass to reduce grass cover to permitted maximum (gC)
real(sp) :: lm_old

real(sp), dimension(npft) :: pft_excess

!-------------------------

present   => vegvars(grid)%present
fpc_grid  => vegvars(grid)%fpc_grid
fpc_ind   => vegvars(grid)%fpc_ind
fpc_inc   => vegvars(grid)%fpc_inc
lm_ind    => vegvars(grid)%lm_ind
rm_ind    => vegvars(grid)%rm_ind
sm_ind    => vegvars(grid)%sm_ind
hm_ind    => vegvars(grid)%hm_ind
nind      => vegvars(grid)%nind
sla       => vegvars(grid)%sla
crownarea => vegvars(grid)%crownarea
lai_ind   => vegvars(grid)%lai_ind

litter_ag_fast => vegvars(grid)%litter_ag_fast
litter_ag_slow => vegvars(grid)%litter_ag_slow
litter_bg      => vegvars(grid)%litter_bg
turnover_ind   => vegvars(grid)%turnover_ind

!-------------------------
! Calculate total woody FPC, FPC increment and grass cover (= crownarea)

where (crownarea > 0.)

  lai_ind  = lm_ind * sla / crownarea
  fpc_ind  = 1. - exp(-0.5 * lai_ind)
  fpc_grid = fpc_ind * nind * crownarea

elsewhere

  lai_ind  = 0.
  fpc_ind  = 0.
  fpc_grid = 0.

end where

ntree           = count(present .and. tree)
fpc_inc_tree    = sum(fpc_inc,  mask = present .and. tree)
fpc_tree_total  = sum(fpc_grid, mask = present .and. tree)

ngrass          = count(present .and. .not. tree)
grasscover      = sum(crownarea, mask = present .and. .not.tree)
fpc_grass_total = sum(fpc_grid,  mask = present .and. .not.tree)

!-------------------------
! Calculate light competition for trees and corresponding fpc_grid reduction

pft_excess = 0.

if (fpc_tree_total > fpc_tree_max) then  ! Reduce tree cover

  excess = fpc_tree_total - fpc_tree_max

  do pft = 1,npft

    if (present(pft) .and. tree(pft) .and. fpc_grid(pft) > 0.) then

      ! This formulation ensures equal competition (precludes total dominance by one PFT)

      pft_excess(pft) = min(fpc_grid(pft), excess * fpc_grid(pft) / fpc_tree_total)

      !original LPJ formulation allows one PFT to become dominant if it has no fpc_inc (so the others are reduced)

      !if (fpc_inc_tree > 0.) then
      !  pft_excess(pft) = min(fpc_grid(pft),excess * (fpc_inc(pft) / fpc_inc_tree))
      !else
      !  pft_excess(pft) = min(fpc_grid(pft),excess / real(ntree))
      !end if

    else

      pft_excess(pft) = 0.

    end if

  end do ! PFT loop

  !-------------------------

  do pft = 1, npft

    if (pft_excess(pft) > 0.) then

      ! Reduce individual density (and thereby gridcell-level biomass) so that total tree FPC reduced to 'fpc_tree_max'

      nind_kill = nind(pft) * pft_excess(pft) / fpc_grid(pft)

      nind(pft) = nind(pft) - nind_kill

      if (lprint .and. grid == gprint) write(0,*) 'NIND kill:', nind_kill

      !Transfer lost biomass to litter

      litter_ag_fast(pft) = litter_ag_fast(pft) + nind_kill * lm_ind(pft)                ! leaf
      litter_ag_slow(pft) = litter_ag_slow(pft) + nind_kill *(sm_ind(pft) + hm_ind(pft)) ! stem
      litter_bg(pft)      = litter_bg(pft)      + nind_kill * rm_ind(pft)                ! root

      !zero out any pft that has been reduced to zero nind

      if (nind(pft) <= 0.) then

        present(pft)  = .false.
        nind(pft)     = 0.
        lm_ind(pft) = 0.
        sm_ind(pft) = 0.
        hm_ind(pft) = 0.
        rm_ind(pft) = 0.

      end if

      !update isotopes
      !ignored for now

    end if ! pft_excess IF condition

  end do ! PFT loop

end if ! fpc_tree_total > fpc_tree_max IF condition


!-------------------------
! Calculate grass competition (inferrior to woody PFTs)

fpc_grass_max = 1. - min(fpc_tree_total, fpc_tree_max)

if (fpc_grass_total > fpc_grass_max) then  ! Reduce grass cover

  excess = fpc_grass_total - fpc_grass_max

  do pft = 1, npft

    if (present(pft) .and. .not.tree(pft) .and. lm_ind(pft) > 0.) then

      lm_old = lm_ind(pft)

      lm_ind(pft) = max(0., -2. * log(1. - (fpc_grid(pft) - excess)) / sla(pft))

      lm_kill = lm_old - lm_ind(pft)
      rm_kill = rm_ind(pft) * lm_kill / lm_old

      rm_ind(pft) = rm_ind(pft) - rm_kill

      ! Transfer lost biomass to litter
      litter_ag_fast(pft) = litter_ag_fast(pft) + lm_kill
      litter_bg(pft)      = litter_bg(pft)      + rm_kill

      !update isotopes
      !ignored for now

    end if

  end do

end if

!-------------------------
! Update LAI and FPC and litter pools after population reduction

where (crownarea > 0.)
  lai_ind  = lm_ind * sla / crownarea
  fpc_ind  = 1. - exp(-0.5 * lai_ind)
  fpc_grid = fpc_ind * nind * crownarea
elsewhere
  lai_ind  = 0.
  fpc_ind  = 0.
  fpc_grid = 0.
end where

! Correct for mathematical overshoot
litter_ag_fast = max(litter_ag_fast, 0.)
litter_ag_slow = max(litter_ag_slow, 0.)
litter_bg      = max(litter_bg, 0.)

! do pft = 1, npft
!   if (fpc_grid(pft) < 0. .or. fpc_grid(pft) > 1.) write(0,*) 'resetting pft ',pft, fpc_grid(pft)
! end do

fpc_grid = max(fpc_grid, 0.)
fpc_grid = min(fpc_grid, 1.)

!-------------------------

end subroutine light


end module lightmod
