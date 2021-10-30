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

use statevarsmod, only : vegvars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: fpc_tree_max = 0.95    ! Maximum tree FPC (Sitch et al., 2003; P. 170)

!-------------------------
! Pointer variables
logical,  pointer :: present           ! PFT present
real(sp), pointer :: fpc_grid          ! Foilage projective cover over grid (fraction)
real(sp), pointer :: fpc_ind           ! Foliage projective cover of individual (fraction)
real(sp), pointer :: fpc_inc           ! Foliage projective cover increment (fraction)
real(sp), pointer :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), pointer :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer :: nind              ! PFT population
real(sp), pointer :: sla               ! Specific leaf area (m2 gC-1)
real(sp), pointer :: crownarea         ! Tree crownarea (m2)
real(sp), pointer :: lai_ind           ! Leaf area index of individual (m2 m-2)
real(sp), pointer :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), pointer :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), pointer :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), pointer :: turnover_ind      ! Total turnover of individual (gC m-2)

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
real(sp) :: pft_excess
real(sp) :: pft_sorted

!-------------------------

present   => vegvars(grid,day)%present
fpc_grid  => vegvars(grid,day)%fpc_grid
fpc_ind   => vegvars(grid,day)%fpc_ind
fpc_inc   => vegvars(grid,day)%fpc_inc
lm_ind    => vegvars(grid,day)%lm_ind
rm_ind    => vegvars(grid,day)%rm_ind
sm_ind    => vegvars(grid,day)%sm_ind
hm_ind    => vegvars(grid,day)%hm_ind
nind      => vegvars(grid,day)%nind
sla       => vegvars(grid,day)%sla
crownarea => vegvars(grid,day)%crownarea
lai_ind   => vegvars(grid,day)%lai_ind

litter_ag_fast => vegvars(grid,day)%litter_ag_fast
litter_ag_slow => vegvars(grid,day)%litter_ag_slow
litter_bg      => vegvars(grid,day)%litter_bg
turnover_ind   => vegvars(grid,day)%turnover_ind

!-------------------------
! Calculate total woody FPC, FPC increment and grass cover (= crownarea)

if (crownarea > 0.) then

  lai_ind  = lm_ind * sla / crownarea
  fpc_ind  = 1. - exp(-0.5 * lai_ind)
  fpc_grid = fpc_ind * nind * crownarea

else

  lai_ind  = 0.
  fpc_ind  = 0.
  fpc_grid = 0.

end if

! ntree           = count(present .and. tree)
! fpc_inc_tree    = sum(fpc_inc,  mask = present .and. tree)
! fpc_tree_total  = sum(fpc_grid, mask = present .and. tree)
!
! ngrass          = count(present .and. .not. tree)
! grasscover      = sum(crownarea, mask = present .and. .not. tree)
! fpc_grass_total = sum(fpc_grid,  mask = present .and. .not. tree)

ntree           = 1
fpc_inc_tree    = fpc_inc
fpc_tree_total  = fpc_grid

ngrass          = 0
grasscover      = 0.
fpc_grass_total = 0.

!-------------------------
! Calculate light competition for trees and corresponding fpc_grid reduction

pft_excess = 0.

if (fpc_tree_total > fpc_tree_max) then  ! Reduce tree cover

  excess = fpc_tree_total - fpc_tree_max

  ! do pft = 1,npft
    if (present .and. fpc_grid > 0.) then

      ! This formulation ensures equal competition (precludes total dominance by one PFT)

      pft_excess = min(fpc_grid, excess *  fpc_grid / fpc_tree_total)

      !original LPJ formulation allows one PFT to become dominant if it has no fpc_inc (so the others are reduced)

      !if (fpc_inc_tree > 0.) then
      !  pft_excess(pft) = min(fpc_grid(pft),excess * (fpc_inc(pft) / fpc_inc_tree))
      !else
      !  pft_excess(pft) = min(fpc_grid(pft),excess / real(ntree))
      !end if

    else

      pft_excess = 0.

    end if

  ! end do

  ! do pft = 1,npft

    if (pft_excess > 0.) then

      ! Reduce individual density (and thereby gridcell-level biomass) so that total tree FPC reduced to 'fpc_tree_max'

      nind_kill = nind * pft_excess / fpc_grid

      if (lprint .and. grid == gprint) write(0,*) 'NIND kill:', nind_kill

      nind = nind - nind_kill

      !Transfer lost biomass to litter

      litter_ag_fast = litter_ag_fast + nind_kill * lm_ind           ! Leaves
      litter_ag_slow = litter_ag_slow + nind_kill *(sm_ind + hm_ind) ! Stems
      litter_bg      = litter_bg      + nind_kill * rm_ind           ! Roots

      !zero out any pft that has been reduced to zero nind

      if (nind <= 0.) then

        present  = .false.
        nind     = 0.
        lm_ind = 0.
        sm_ind = 0.
        hm_ind = 0.
        rm_ind = 0.

      end if

      !update isotopes
      !ignored for now

    end if

  ! end do

end if

! #############################
! GRASS COMPETITION HERE
! #############################

!-------------------------
! Update LAI and FPC and litter pools after population reduction

if (crownarea > 0.) then
  lai_ind  = lm_ind * sla / crownarea
  fpc_ind  = 1. - exp(-0.5 * lai_ind)
  fpc_grid = fpc_ind * nind * crownarea
else
  lai_ind  = 0.
  fpc_ind  = 0.
  fpc_grid = 0.
end if

! Correct for mathematical overshoot
litter_ag_fast = max(litter_ag_fast,0.)
litter_ag_slow = max(litter_ag_slow,0.)
litter_bg      = max(litter_bg,0.)

fpc_grid = max(fpc_grid,0.)
fpc_grid = min(fpc_grid,1.)

!-------------------------

end subroutine light


end module lightmod
