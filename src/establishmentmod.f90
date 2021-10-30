module establishmentmod

! Module to calculate vegation establishment
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

real(sp) :: lm_sapl
real(sp) :: hm_sapl
real(sp) :: sm_sapl
real(sp) :: rm_sapl
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

use statevarsmod, only : vegvars

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

real(sp), parameter :: pi    = 3.1415926535
real(sp), parameter :: allom1 = 100
real(sp), parameter :: allom2 = 40
real(sp), parameter :: allom3 = 0.5
real(sp), parameter :: allom4 = 0.3
real(sp), parameter :: leaf_longevity = 1.0
real(sp), parameter :: lai_sapl = 1.5
real(sp), parameter :: reinickerp = 1.6
real(sp), parameter :: latosa     = 8.e3
real(sp), parameter :: wooddens   = 2.e5

real(sp), pointer :: sla

real(sp) :: x
real(sp) :: lmtorm

sla => vegvars(grid,day)%sla

x = 3.

sla = 2e-4 * exp(6.15) / ((12. * leaf_longevity) ** 0.46)

lm_sapl = (lai_sapl*allom1*x**reinickerp*(4.0*sla/pi/latosa)**(reinickerp*0.5) / sla)**(1.0-1.0/reinickerp)

stemdiam_sapl = x*(4.0*lm_sapl*sla/pi/latosa)**0.5

height_sapl = allom2 * stemdiam_sapl**allom3   !Eqn 16

sm_sapl = wooddens * height_sapl * lm_sapl * sla / latosa

hm_sapl = (x-1.0)*sm_sapl  !Eqn 22

lmtorm = 1.0

rm_sapl = (1.0 / lmtorm) * lm_sapl  !From Eqn 23

! print*, lm_sapl, rm_sapl, hm_sapl, sm_sapl, height_sapl, stemdiam_sapl, sla

end subroutine sapling

!---------------------------------------------------------------------

subroutine establishment(year,grid,day)

use statevarsmod, only : ndyear,dayvars,soilvars,vegvars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: pi    = 3.1415926535
real(sp), parameter :: allom1 = 100.
real(sp), parameter :: allom2 = 40.
real(sp), parameter :: allom3 = 0.5
real(sp), parameter :: allom4 = 0.3
real(sp), parameter :: reinickerp = 1.6
real(sp), parameter :: latosa     = 8.e3
real(sp), parameter :: wooddens   = 2.e5
real(sp), parameter :: aprec_min_estab  = 100.        !minimum annual precipitation for establishment (mm)
real(sp), parameter :: estab_max        = 0.24        !maximum sapling establishment rate (indiv/m2)
real(sp), parameter :: nind_min         = 1.e-10      !minimum individual density for persistence of PFT (indiv/m2)
real(sp), parameter :: eps              = 1.e-6       !epsilon parameter (?)

real(sp), parameter :: crownarea_max = 15.

!-------------------------
! Pointer variables
logical, pointer :: present
logical, pointer :: estab
logical, pointer :: survive
real(sp), pointer :: dwscal
real(sp), pointer :: fpc_grid
real(sp), pointer :: fpc_ind
real(sp), pointer :: fpc_inc
real(sp), pointer :: nind
real(sp), pointer :: sla
real(sp), pointer :: stemdiam
real(sp), pointer :: height
real(sp), pointer :: crownarea
real(sp), pointer :: lai_ind
real(sp), pointer :: hm_ind
real(sp), pointer :: lm_ind
real(sp), pointer :: sm_ind
real(sp), pointer :: rm_ind

!-------------------------
! Local variables
real(sp) :: aprec
real(sp) :: estab_grid
real(sp) :: estab_mass
real(sp) :: estab_rate
real(sp) :: estab_pft
real(sp) :: fpc_total
real(sp) :: fpc_tree_total
real(sp) :: nind0
real(sp) :: sm_ind_tmp
integer(i4) :: npft_estab
integer(i4) :: ngrass

!-------------------------

present => vegvars(grid,day)%present
estab => vegvars(grid,day)%estab
survive => vegvars(grid,day)%survive
dwscal => vegvars(grid,day)%dwscal
fpc_grid => vegvars(grid,day)%fpc_grid
fpc_ind => vegvars(grid,day)%fpc_ind
fpc_inc => vegvars(grid,day)%fpc_inc
nind => vegvars(grid,day)%nind
sla => vegvars(grid,day)%sla
stemdiam => vegvars(grid,day)%stemdiam
height => vegvars(grid,day)%height
crownarea => vegvars(grid,day)%crownarea
lai_ind => vegvars(grid,day)%lai_ind
hm_ind => vegvars(grid,day)%hm_ind
lm_ind => vegvars(grid,day)%lm_ind
sm_ind => vegvars(grid,day)%sm_ind
rm_ind => vegvars(grid,day)%rm_ind

!-------------------------

if (year == 1) then
  present = .false.
  survive = .true.
  estab = .true.
  dwscal = 1.0
  fpc_grid = 0.0
  fpc_ind = 0.0
  fpc_inc = 0.0
  nind = 0.0
  height = 2.
  crownarea = 0.0
  lai_ind = 0.0
  hm_ind = 0.0
  lm_ind = 0.0
  sm_ind = 0.0
  rm_ind = 0.0
end if

!-------------------------
! Introduce PFT if conditions suitable for establishment

aprec = sum(dayvars(grid,:)%prec)

if (.not.present .and. survive .and. estab .and. aprec > aprec_min_estab) then

  present = .true.
  nind = 0.0

end if

!-------------------------
! Calculate total woody FPC and ability to establish

ngrass = 0
if(present) npft_estab = 1
fpc_total = fpc_grid
fpc_tree_total = fpc_grid

if (aprec >= aprec_min_estab .and. npft_estab > 0) then

  ! Calculate establishment rate over available space, per tree PFT
  ! Maximum establishment rate reduced by shading as tree FPC approches 1
  ! NOTE : Total establihsment rate partitioned equqally among regenerating woody PFTs
  estab_rate = estab_max * (1. - exp(-5. * (1. - fpc_tree_total))) / real(npft_estab)

  estab_grid = max(estab_rate * (1. - fpc_tree_total),0.)

  ! TEMPORARY condition (Leo Oct 2021)
  if (fpc_grid >= 1.) estab_grid = 0.

else  ! Unsuitable climate for establishment

  estab_grid = 0.

end if

if (lprint .and. grid==gprint)  print *, fpc_grid, fpc_tree_total, estab_rate, estab_grid

!-------------------------
! Establish new sapling carbon mass to current population

if (estab_grid > 0.0) then

  nind0 = nind    ! Save copy of current population
  nind = nind0 + estab_grid

  sm_ind_tmp = (sm_ind * nind0 + sm_sapl * estab_grid) / nind ! Sapwood mass
  hm_ind = (hm_ind * nind0 + hm_sapl * estab_grid) / nind ! Heartwood mass
  lm_ind = (lm_ind * nind0 + lm_sapl * estab_grid) / nind ! Leaf mass
  rm_ind = (rm_ind * nind0 + rm_sapl * estab_grid) / nind ! Root mass

  !Accumulate biomass increment due to sapling establishment

  estab_mass = lm_sapl + sm_sapl + hm_sapl + rm_sapl

  estab_pft = estab_mass * estab_grid

  ! if (estab_mass * estab_grid > eps) acflux_estab(1) = acflux_estab(1) + estab_mass * estab_grid

  stemdiam = (4.0 * (sm_ind_tmp + hm_ind) / wooddens / pi / allom2) ** (1.0 / (2.0 + allom3)) !Eqn 9

  height = allom2 * (stemdiam ** allom3)                          !Eqn C

  crownarea = min(crownarea_max, allom1 * (stemdiam ** reinickerp)) !Eqn D

  !Recalculate sapwood mass, transferring excess sapwood to heartwood compartment, if necessary to satisfy Eqn A

  sm_ind = lm_ind * height * wooddens * sla / latosa

  hm_ind = max(hm_ind + (sm_ind_tmp - sm_ind),0.)

end if

!-------------------------
! Update LAI and FPC

if (present) then

  if (crownarea > 0.0) then
    lai_ind = (lm_ind * sla) / crownarea
  else
    lai_ind = 0.0
  end if

  fpc_ind = (1. - exp(-0.5 * lai_ind))
  fpc_grid = crownarea * nind * fpc_ind

end if

! if (lprint .and. grid==gprint) print *, lm_ind, rm_ind, sm_ind, hm_ind, nind, &
!                                         height, stemdiam, crownarea, lai_ind, fpc_ind, fpc_grid

end subroutine establishment

!---------------------------------------------------------------------

end module establishmentmod
