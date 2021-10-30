module establishmentmod

! Module to calculate vegation establishment
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

!-------------------------
! Module parameters
real(sp), parameter :: pi             = 3.1415926535
real(sp), parameter :: allom1         = 100           ! Allometric constant 1 (Sitch et al., 2003; Table 3)
real(sp), parameter :: allom2         = 40            ! Allometric constant 2 (Sitch et al., 2003; Table 3)
real(sp), parameter :: allom3         = 0.5           ! Allometric constant 3 (Sitch et al., 2003; Table 3)
real(sp), parameter :: allom4         = 0.3           ! Allometric constant 4
real(sp), parameter :: leaf_longevity = 1.0           ! a_leaf (Sitch et al., 2003; Table 1)
real(sp), parameter :: lai_sapl       = 1.5           ! Sapling leaf area index (Sitch et al., 2003; P.171)
real(sp), parameter :: reinickerp     = 1.6           ! Reinecke's constant (Sitch et al., 2003; Eq. 4)
real(sp), parameter :: latosa         = 8.e3          ! Individual leaf area to sapwood cross sectional area ratio (Sitch et al., 2003; Table 3)
real(sp), parameter :: wooddens       = 2.e5          ! Sapwood wood density (kg m-3)

!-------------------------
! Module variables for sapling allometry
real(sp) :: lm_sapl
real(sp) :: rm_sapl
real(sp) :: sm_sapl
real(sp) :: hm_sapl
real(sp) :: stemdiam_sapl
real(sp) :: height_sapl

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine bioclim(year,grid,day)

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day


end subroutine bioclim

!---------------------------------------------------------------------

subroutine sapling(year,grid,day)

! Subroutine to initiate sapling parameters
! Saplings at establishment are assumed constant allometry throughout model

use statevarsmod, only : vegvars

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: lmtorm = 1.   ! Leaf-to-root mass ratio (Sitch et al., 2003; Table 3)
real(sp), parameter :: x      = 3.   ! Sapling (sapwood + heartwood) / sapwood mass ratio

real(sp), pointer :: sla

!-------------------------

sla => vegvars(grid,day)%sla

!-------------------------

! Calculate specific leaf area (m2 gC-1)
! Sitch et al. (2003) Eq. 6
sla = 2e-4 * exp(6.15) / ((12. * leaf_longevity) ** 0.46)

! Calculate leaf mass of sapling individual (gC)
! Modified from Sitch et al. (2003) Eq. 5
lm_sapl = (lai_sapl * allom1 * x ** reinickerp * &
          (4.0 * sla / pi / latosa) ** (reinickerp * 0.5) / sla)**(1. - 1. / reinickerp)

! Calculate stem diameter of sapling individual (m)
stemdiam_sapl =  x * (4. * lm_sapl * sla / pi / latosa) ** 0.5

! Calculate height of sapling individual (m)
! Sitch et al. (2003) Eq. 3
height_sapl = allom2 * stemdiam_sapl ** allom3

! Calculate sapwood mass of sapling individual (gC)
sm_sapl = wooddens * height_sapl * lm_sapl * sla / latosa

! Calculate heartwood mass of sapling individual (gC)
hm_sapl = (x - 1.0) * sm_sapl

! Calculate root mass of sapling individual (gC)
rm_sapl = (1.0 / lmtorm) * lm_sapl


end subroutine sapling

!---------------------------------------------------------------------

subroutine establishment(year,grid,day)

! Subroutine to establish sapling with grid have available space

use statevarsmod, only : ndyear,dayvars,soilvars,vegvars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: aprec_min_estab  = 100.      ! Minimum annual precipitation for establishment (mm)
real(sp), parameter :: estab_max        = 0.24      ! Maximum sapling establishment rate (indiv/m2) (Sitch et al., 2003; P.171)
real(sp), parameter :: nind_min         = 1.e-10    ! Minimum individual density for persistence of PFT (indiv/m2)
real(sp), parameter :: eps              = 1.e-6     ! Epsilon parameter (?)
real(sp), parameter :: crownarea_max    = 15.       ! Maximum crownarea for individual (m2)

!-------------------------
! Pointer variables
logical,  pointer :: present           ! PFT present
logical,  pointer :: estab             ! PFT establishment
logical,  pointer :: survive           ! PFT survival
real(sp), pointer :: dwscal            ! Daily water stress factor (supply/demand ratio)
real(sp), pointer :: fpc_grid          ! Foilage projective cover over grid (fraction)
real(sp), pointer :: fpc_ind           ! Foliage projective cover of individual (fraction)
real(sp), pointer :: fpc_inc           ! Foliage projective cover increment (fraction)
real(sp), pointer :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp), pointer :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp), pointer :: nind              ! PFT population
real(sp), pointer :: sla               ! Specific leaf area (m2 gC-1)
real(sp), pointer :: stemdiam          ! Tree stem diameter (m)
real(sp), pointer :: height            ! Tree height (m)
real(sp), pointer :: crownarea         ! Tree crownarea (m2)
real(sp), pointer :: lai_ind           ! Leaf area index of individual (m2 m-2)

!-------------------------
! Local variables
real(sp)    :: aprec            ! Annual total precipitation (mm yr-1)
real(sp)    :: estab_grid       ! Grid-level establishment rate (indiv/m2)
real(sp)    :: estab_rate       ! Sapling establishment rate over area available for establishment (indiv/m2)
real(sp)    :: estab_mass       ! Mass of establishing PFT (gC ind-1)
real(sp)    :: estab_pft
real(sp)    :: fpc_total        ! Total grid fpc (fraction)
real(sp)    :: fpc_tree_total   ! Total grid FPC for tree PFT (fraction)
real(sp)    :: nind0            ! Number of individual before establishment (indiv)
real(sp)    :: sm_ind_tmp       ! Temporary sapwood variable (gC m-2)
integer(i4) :: npft_estab       ! Number of regenerating tree PFTs
integer(i4) :: ngrass           ! Number of grass PFT present

!-------------------------

present   => vegvars(grid,day)%present
estab     => vegvars(grid,day)%estab
survive   => vegvars(grid,day)%survive

fpc_grid  => vegvars(grid,day)%fpc_grid
fpc_ind   => vegvars(grid,day)%fpc_ind
fpc_inc   => vegvars(grid,day)%fpc_inc
dwscal    => vegvars(grid,day)%dwscal

lm_ind    => vegvars(grid,day)%lm_ind
rm_ind    => vegvars(grid,day)%rm_ind
sm_ind    => vegvars(grid,day)%sm_ind
hm_ind    => vegvars(grid,day)%hm_ind
nind      => vegvars(grid,day)%nind
sla       => vegvars(grid,day)%sla
stemdiam  => vegvars(grid,day)%stemdiam
height    => vegvars(grid,day)%height
crownarea => vegvars(grid,day)%crownarea
lai_ind   => vegvars(grid,day)%lai_ind

!-------------------------
! Initialize values for vegvars

if (year == 1) then
  present   = .false.
  survive   = .true.
  estab     = .true.
  dwscal    = 1.
  fpc_grid  = 0.
  fpc_ind   = 0.
  fpc_inc   = 0.
  nind      = 0.
  height    = 2.
  crownarea = 0.
  lai_ind   = 0.
  lm_ind    = 0.
  rm_ind    = 0.
  sm_ind    = 0.
  hm_ind    = 0.
end if

!-------------------------
! Introduce PFT if conditions suitable for establishment

aprec = sum(dayvars(grid,:)%prec)

if (.not.present .and. survive .and. estab .and. aprec > aprec_min_estab) then

  present = .true.
  nind    = 0.0

end if

!-------------------------
! Calculate total woody FPC and ability to establish

ngrass = 0                    ! No grass PFT at the moment (Leo Lai Oct 2021)
if(present) npft_estab = 1    ! Currently only consider one generic "PFT" (Leo Lai Oct 2021)
fpc_total      = fpc_grid
fpc_tree_total = fpc_grid

if (aprec >= aprec_min_estab .and. npft_estab > 0) then

  ! Calculate establishment rate over available space, per tree PFT
  ! Maximum establishment rate reduced by shading as tree FPC approches 1
  ! NOTE : Total establihsment rate partitioned equqally among regenerating woody PFTs
  estab_rate = estab_max * (1. - exp(-5. * (1. - fpc_tree_total))) / real(npft_estab)

  estab_grid = max(estab_rate * (1. - fpc_tree_total),0.)

else  ! Unsuitable climate for establishment

  estab_grid = 0.

end if

! if (lprint .and. grid==gprint)  print *, fpc_grid, fpc_tree_total, estab_rate, estab_grid

!-------------------------
! Establish new sapling carbon mass to current population

if (estab_grid > 0.0) then

  ! Save copy of current population and accumulate new sapling population
  nind0 = nind
  nind = nind0 + estab_grid

  ! Calculate new carbon mass pools after accumulating new population
  lm_ind     = (lm_ind * nind0 + lm_sapl * estab_grid) / nind ! Leaf mass
  rm_ind     = (rm_ind * nind0 + rm_sapl * estab_grid) / nind ! Root mass
  sm_ind_tmp = (sm_ind * nind0 + sm_sapl * estab_grid) / nind ! Sapwood mass
  hm_ind     = (hm_ind * nind0 + hm_sapl * estab_grid) / nind ! Heartwood mass

  !Accumulate biomass increment due to sapling establishment
  estab_mass = lm_sapl + sm_sapl + hm_sapl + rm_sapl

  estab_pft = estab_mass * estab_grid

  ! if (estab_mass * estab_grid > eps) acflux_estab(1) = acflux_estab(1) + estab_mass * estab_grid

  ! Calculate new allometric variables from new carbon pools
  stemdiam = (4.0 * (sm_ind_tmp + hm_ind) / wooddens / pi / allom2) ** (1.0 / (2.0 + allom3))

  height = allom2 * (stemdiam ** allom3)

  crownarea = min(crownarea_max, allom1 * (stemdiam ** reinickerp))

  ! Recalculate sapwood mass, transferring excess sapwood to heartwood compartment, if necessary to satisfy Eqn A
  sm_ind = lm_ind * height * wooddens * sla / latosa

  hm_ind = max(hm_ind + (sm_ind_tmp - sm_ind), 0.)

end if

!-------------------------
! Update LAI and FPC
if (present) then

  if (crownarea > 0.0) then
    lai_ind = (lm_ind * sla) / crownarea
  else
    lai_ind = 0.0
  end if

  !------

  fpc_ind  = (1. - exp(-0.5 * lai_ind))
  fpc_grid = crownarea * nind * fpc_ind

end if

!-------------------------

! if (lprint .and. grid==gprint) print *, lm_ind, rm_ind, sm_ind, hm_ind, nind, &
!                                         height, stemdiam, crownarea, lai_ind, fpc_ind, fpc_grid

end subroutine establishment

!---------------------------------------------------------------------

end module establishmentmod
