module establishmentmod

! Module to calculate vegation establishment
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp
use pftparmod,     only : npft

implicit none

!-------------------------
! Module parameters
real(sp), parameter :: pi             = 3.1415926535
real(sp), parameter :: allom1         = 100.          ! Allometric constant 1 (Sitch et al., 2003; Table 3)
real(sp), parameter :: allom2         = 40.           ! Allometric constant 2 (Sitch et al., 2003; Table 3)
real(sp), parameter :: allom3         = 0.5           ! Allometric constant 3 (Sitch et al., 2003; Table 3)
real(sp), parameter :: allom4         = 0.3           ! Allometric constant 4
real(sp), parameter :: reinickerp     = 1.6           ! Reinecke's constant (Sitch et al., 2003; Eq. 4)
! real(sp), parameter :: latosa         = 8.e3          ! Individual leaf area to sapwood cross sectional area ratio (Sitch et al., 2003; Table 3)
! real(sp), parameter :: wooddens       = 2.e5          ! Sapwood wood density (kg m-3)

!-------------------------
! Module variables for sapling allometry
real(sp), dimension(npft) :: lm_sapl
real(sp), dimension(npft) :: rm_sapl
real(sp), dimension(npft) :: sm_sapl
real(sp), dimension(npft) :: hm_sapl
real(sp), dimension(npft) :: crownarea_sapl
real(sp), dimension(npft) :: stemdiam_sapl
real(sp), dimension(npft) :: height_sapl

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine bioclim(year,grid,spinup,ndyear,mtemp,dtemp,estab,survive)

use pftparmod,    only : npft,pftpar
use statevarsmod, only : invars

integer(i4),               intent(in)    :: year
integer(i4),               intent(in)    :: grid
logical,                   intent(in)    :: spinup
integer(i4),               intent(in)    :: ndyear
real(sp),    dimension(:), intent(in)    :: mtemp             ! Annual 12month monthly temperature series (degC)
real(sp),    dimension(:), intent(in)    :: dtemp             ! Annual daily tmean series (degC)
logical,     dimension(:), intent(inout) :: estab             ! PFT establishment
logical,     dimension(:), intent(inout) :: survive           ! PFT survival

real(sp), parameter :: gddbase = 5.

real(sp), allocatable, dimension(:)  :: mtemp_20

real(sp) :: mtemp_min20
real(sp) :: mtemp_max
real(sp) :: gdd

real(sp) :: tcmin   ! PFT-specific minimum coldest-month temperature
real(sp) :: tcmax   ! PFT-specific maximum coldest-month temperature
real(sp) :: gddmin  ! PFT-specific minimum GDD
real(sp) :: twmax   ! PFT-specific upper limit of warmest-month temperature

integer :: pft
integer :: s
integer :: e
integer :: ny
integer :: i

!-------------------------
! Initialize values for vegvars

if (year == 1 .and. spinup) then
  survive   = .false.
  estab     = .false.
end if

!-------------------------
! Calculate average monthly minimum temperature for the last 20 years

if (year < 20) then

  s = 13
  e = 12 * (year + 1)

  allocate(mtemp_20(12*year))

  mtemp_20 = invars(grid)%tmp(s:e)
  ny = year

else

  e = 12 * (year + 1)
  s = e - (20 * 12) + 1

  allocate(mtemp_20(240))

  mtemp_20 = invars(grid)%tmp(s:e)
  ny = 20

end if

!------

s = 1
e = 12
mtemp_min20 = 0.

do i = 1, ny

   mtemp_min20 = mtemp_min20 + minval(mtemp_20(s:e))

   s = s + 12
   e = e + 12

end do

mtemp_min20 = mtemp_min20 / ny

!-------------------------
! Calculate maximum monthly temperature of the year and GDD5

mtemp_max = maxval(mtemp)

gdd = sum(dtemp-gddbase, mask = dtemp > gddbase)

!-------------------------

do pft = 1, npft

  !------

  tcmin  = pftpar%tcmin(pft)
  tcmax  = pftpar%tcmax(pft)
  gddmin = pftpar%gddmin(pft)
  twmax  = pftpar%twmax(pft)

  !------

  if (mtemp_min20 >= tcmin) then

    survive(pft) = .true.

    if (gdd >= gddmin .and. mtemp_min20 <= tcmax .and. mtemp_max <= twmax) then

      estab(pft) = .true.

    else

      estab(pft) = .false.

    end if

  else

    survive(pft) = .false.

  end if

end do

deallocate(mtemp_20)


end subroutine bioclim

!---------------------------------------------------------------------

subroutine sapling(sla)

! Subroutine to initiate sapling parameters
! Saplings at establishment are assumed constant allometry throughout model

use parametersmod, only : i4,sp
use pftparmod,     only : npft,pftpar,tree

real(sp), dimension(:), intent(out) :: sla

! Local variables
integer(i4) :: pft
real(sp)    :: leaf_longevity
real(sp)    :: lai_sapl
real(sp)    :: x
real(sp)    :: lmtorm
real(sp)    :: latosa
real(sp)    :: wooddens

!-------------------------

do pft = 1, npft

  !------

  leaf_longevity = pftpar%longevity(pft)
  lai_sapl       = pftpar%lai_sapl(pft)
  x              = pftpar%x_sapl(pft)
  lmtorm         = pftpar%lmtorm(pft)
  latosa         = pftpar%latosa(pft)
  wooddens       = pftpar%wooddens(pft)

  !------

  ! Calculate specific leaf area (m2 gC-1)
  ! Sitch et al. (2003) Eq. 6
  sla(pft) = 2.e-4 * exp(6.15) / ((12. * leaf_longevity) ** 0.46)

  if (tree(pft)) then

    ! Calculate leaf mass of sapling individual (gC)
    ! Modified from Sitch et al. (2003) Eq. 5
    lm_sapl(pft) = (lai_sapl * allom1 * x ** reinickerp * &
                  (4.0 * sla(pft) / pi / latosa) ** (reinickerp * 0.5) / sla(pft))**(1. - 1. / reinickerp)

    ! Calculate the crownarea for sapling individual (m2)
    crownarea_sapl(pft) = (lm_sapl(pft) * sla(pft)) / lai_sapl

    ! Calculate stem diameter of sapling individual (m)
    stemdiam_sapl(pft) =  x * (4. * lm_sapl(pft) * sla(pft) / pi / latosa) ** 0.5

    ! Calculate height of sapling individual (m)
    ! Sitch et al. (2003) Eq. 3
    height_sapl(pft) = allom2 * stemdiam_sapl(pft) ** allom3

    ! Calculate sapwood mass of sapling individual (gC)
    sm_sapl(pft) = wooddens * height_sapl(pft) * lm_sapl(pft) * sla(pft) / latosa

    ! Calculate heartwood mass of sapling individual (gC)
    hm_sapl(pft) = (x - 1.0) * sm_sapl(pft)

    ! Calculate root mass of sapling individual (gC)
    rm_sapl(pft) = (1.0 / lmtorm) * lm_sapl(pft)

  else ! Grass

    ! Leaf mass
    lm_sapl(pft) = lai_sapl / sla(pft)

    ! Calculate initial root mass
    rm_sapl(pft) = (1.0 / lmtorm) * lm_sapl(pft)

    ! No saphood or heartwood mass for grass PFTs
    ! sm_sapl = 0.
    ! hm_sapl = 0.

  end if

end do


end subroutine sapling

!---------------------------------------------------------------------

subroutine establishment(year,spinup,dprec,present,estab,survive,fpc_grid,fpc_ind,fpc_inc,lm_ind,rm_ind,sm_ind,hm_ind, &
                         nind,sla,stemdiam,height,crownarea,lai_ind,dwscal,leafon,leafondays,leafoffdays,       &
                         litter_ag_fast,litter_ag_slow,litter_bg)

! Subroutine to establish sapling with grid have available space

use utilitiesmod, only : getmonth
use pftparmod,    only : npft,pftpar,tree

! Arguments
integer(i4),               intent(in)    :: year
logical,                   intent(in)    :: spinup
real(sp),    dimension(:), intent(in)    :: dprec             ! Daily precipitation series (mm)
logical,     dimension(:), intent(inout) :: present           ! PFT present
logical,     dimension(:), intent(inout) :: estab             ! PFT establishment
logical,     dimension(:), intent(inout) :: survive           ! PFT survival
real(sp),    dimension(:), intent(inout) :: fpc_grid          ! Foilage projective cover over grid (fraction)
real(sp),    dimension(:), intent(inout) :: fpc_ind           ! Foliage projective cover of individual (fraction)
real(sp),    dimension(:), intent(inout) :: fpc_inc           ! Foliage projective cover increment (fraction)
real(sp),    dimension(:), intent(inout) :: lm_ind            ! Leaf carbon mass of individual (gC m-2)
real(sp),    dimension(:), intent(inout) :: rm_ind            ! Root carbon mass of individual (gC m-2)
real(sp),    dimension(:), intent(inout) :: sm_ind            ! Sapwood carbon mass of individual (gC m-2)
real(sp),    dimension(:), intent(inout) :: hm_ind            ! Heartwood carbon mass of individual (gC m-2)
real(sp),    dimension(:), intent(inout) :: nind              ! PFT population
real(sp),    dimension(:), intent(inout) :: sla               ! Specific leaf area (m2 gC-1)
real(sp),    dimension(:), intent(inout) :: stemdiam          ! Tree stem diameter (m)
real(sp),    dimension(:), intent(inout) :: crownarea         ! Tree crownarea (m2)
real(sp),    dimension(:), intent(inout) :: height            ! Tree height (m)
real(sp),    dimension(:), intent(inout) :: lai_ind           ! Leaf area index of individual (m2 m-2)
real(sp),    dimension(:), intent(inout) :: dwscal            ! Daily water stress factor (supply/demand ratio)
logical,     dimension(:), intent(inout) :: leafon
real(sp),    dimension(:), intent(inout) :: leafondays
real(sp),    dimension(:), intent(inout) :: leafoffdays
real(sp),    dimension(:), intent(inout) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp),    dimension(:), intent(inout) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp),    dimension(:), intent(inout) :: litter_bg         ! Below ground litter pool (gC m-2)

!-------------------------
! Parameters
real(sp), parameter :: aprec_min_estab  = 100.0     ! 100.0 mm - Minimum annual precipitation for establishment (mm)
real(sp), parameter :: estab_max        = 0.15      ! Maximum sapling establishment rate (indiv/m2) (Sitch et al., 2003; P.171)
real(sp), parameter :: nind_min         = 1.e-10    ! Minimum individual density for persistence of PFT (indiv/m2)
real(sp), parameter :: eps              = 1.e-6     ! Epsilon parameter (?)

!-------------------------
! Local variables
real(sp)    :: mprec
real(sp)    :: aprec            ! Annual total precipitation (mm yr-1)
real(sp)    :: estab_grid       ! Grid-level establishment rate (indiv/m2)
real(sp)    :: estab_rate       ! Sapling establishment rate over area available for establishment (indiv/m2)
real(sp)    :: estab_mass       ! Mass of establishing PFT (gC ind-1)
real(sp)    :: estab_pft
real(sp)    :: fpc_total        ! Total grid fpc (fraction)
real(sp)    :: fpc_tree_total   ! Total grid FPC for tree PFT (fraction)
real(sp)    :: bare             ! Bare ground not occupied by tree PFT (fraction)
real(sp)    :: nind0            ! Number of individual before establishment (indiv)
real(sp)    :: sm_ind_tmp       ! Temporary sapwood variable (gC m-2)
integer(i4) :: npft_estab       ! Number of regenerating tree PFTs
integer(i4) :: ngrass           ! Number of grass PFT present

real(sp) :: latosa
real(sp) :: wooddens
real(sp) :: crownarea_max
real(sp) :: acflux_estab
real(sp) :: sap_xsa
integer :: pft

integer :: m
integer :: startday
integer :: endday

!-------------------------
! Initialize values for vegvars

if (year == 1 .and. spinup) then
  present   = .false.
  dwscal    = 1.
  fpc_grid  = 0.
  fpc_ind   = 0.
  fpc_inc   = 0.
  nind      = 0.
  height    = 0.
  crownarea = 0.
  lai_ind   = 0.
  lm_ind    = 0.
  rm_ind    = 0.
  sm_ind    = 0.
  hm_ind    = 0.
end if

!-------------------------
! Introduce PFT if conditions suitable for establishment and intialize vegetation statevars

! call getmonth(day,ndyear,m,startday,endday)

aprec = sum(dprec)
! mprec = genvars%tmp(4+m)

do pft = 1, npft
! do pft = 8, npft

  if (present(pft) .and. (.not.survive(pft) .or. nind(pft) < nind_min)) then

    ! Kill PFT
    present(pft) = .false.

    litter_ag_fast(pft) = litter_ag_fast(pft) + nind(pft) * lm_ind(pft)
    litter_ag_slow(pft) = litter_ag_slow(pft) + nind(pft) *(sm_ind(pft) + hm_ind(pft))
    litter_bg(pft)      = litter_bg(pft)      + nind(pft) * rm_ind(pft)

    nind(pft) = 0.

    lm_ind(pft) = 0.
    sm_ind(pft) = 0.
    hm_ind(pft) = 0.
    rm_ind(pft) = 0.

  else if (.not.present(pft) .and. survive(pft) .and. estab(pft) .and. aprec > aprec_min_estab) then
  ! else if (.not.present(pft) .and. survive(pft) .and. estab(pft) .and. mprec > aprec_min_estab/12.) then

    present(pft) = .true.

    if (tree(pft)) then
      nind(pft) = 0.
    else
      nind(pft) = 1.    ! Each grass PFT = 1 "individual"
    end if

    fpc_grid(pft) = 0.

    if (.not.tree(pft)) crownarea(pft) = 1.

    leafon(pft)      =.true.
    leafondays(pft)  = 0
    leafoffdays(pft) = 0

  end if

end do

!-------------------------
! Calculate total woody FPC and ability to establish

ngrass         = count(present .and. .not.tree)
npft_estab     = count(present .and. estab .and. tree)
fpc_total      = sum(fpc_grid)
fpc_tree_total = sum(fpc_grid, mask = tree)

acflux_estab = 0.

if (aprec >= aprec_min_estab .and. npft_estab > 0) then
! if (mprec >= aprec_min_estab/12. .and. npft_estab > 0) then

  ! Calculate establishment rate over available space, per tree PFT
  ! Maximum establishment rate reduced by shading as tree FPC approches 1
  ! NOTE : Total establihsment rate partitioned equqally among regenerating woody PFTs
  estab_rate = estab_max * (1. - exp(-5. * (1. - fpc_tree_total))) / real(npft_estab)

  estab_grid = max(estab_rate * (1. - fpc_tree_total),0.)

else  ! Unsuitable climate for establishment

  estab_grid = 0.

end if

!-------------------------
! Establish new sapling carbon mass to current population

do pft = 1, npft
! do pft = 8, npft

  if (present(pft) .and. estab(pft)) then

    if (tree(pft)) then

      if (estab_grid > 0.0) then

        crownarea_max = pftpar%maxcrowna(pft)
        latosa        = pftpar%latosa(pft)
        wooddens      = pftpar%wooddens(pft)

        ! Save copy of current population and accumulate new sapling population
        nind0 = nind(pft)
        nind(pft) = nind0 + estab_grid

        ! Calculate new carbon mass pools after accumulating new population
        lm_ind(pft) = (lm_ind(pft) * nind0 + lm_sapl(pft) * estab_grid) / nind(pft) ! Leaf mass
        rm_ind(pft) = (rm_ind(pft) * nind0 + rm_sapl(pft) * estab_grid) / nind(pft) ! Root mass
        sm_ind_tmp  = (sm_ind(pft) * nind0 + sm_sapl(pft) * estab_grid) / nind(pft) ! Sapwood mass
        hm_ind(pft) = (hm_ind(pft) * nind0 + hm_sapl(pft) * estab_grid) / nind(pft) ! Heartwood mass

        !Accumulate biomass increment due to sapling establishment
        estab_mass = lm_sapl(pft) + sm_sapl(pft) + hm_sapl(pft) + rm_sapl(pft)

        ! estab_pft(pft) = estab_mass * estab_grid

        ! if (estab_mass * estab_grid > eps) acflux_estab(1) = acflux_estab(1) + estab_mass * estab_grid

        ! Calculate new allometric variables from new carbon pools
        stemdiam(pft) = (4. * (sm_ind_tmp + hm_ind(pft)) / wooddens / pi / allom2) ** (1. / (2. + allom3))

        height(pft) = allom2 * (stemdiam(pft) ** allom3)

        crownarea(pft) = min(crownarea_max, allom1 * (stemdiam(pft) ** reinickerp))

        ! TESTING DIFFERENT SET OF EQUATIONS
        ! sap_xsa = lm_ind(pft) * sla(pft) / latosa
        !
        ! height(pft) = sm_ind_tmp / sap_xsa / wooddens
        !
        ! stemdiam(pft) = (height(pft) / allom2) ** (1. / allom3)

        ! Recalculate sapwood mass, transferring excess sapwood to heartwood compartment, if necessary to satisfy Eqn A
        sm_ind(pft) = lm_ind(pft) * height(pft) * wooddens * sla(pft) / latosa

        hm_ind(pft) = max(hm_ind(pft) + (sm_ind_tmp - sm_ind(pft)), 0.)

      end if ! estab_grid IF condition

    else ! Grass establishment

      ! Grasses can establish in non-vegetated areas
      bare = (1. - fpc_total) / real(ngrass)

      lm_ind(pft) = lm_ind(pft) + bare * lm_sapl(pft)
      rm_ind(pft) = rm_ind(pft) + bare * rm_sapl(pft)

      !Accumulate biomass increment due to grass establishment
      estab_mass = bare * (lm_sapl(pft) + rm_sapl(pft))

      ! estab_pft(pft) = estab_mass * crownarea(pft)

      ! if (estab_mass * crownarea(pft) > eps) acflux_estab(1) = acflux_estab(1) + estab_mass * crownarea(pft)

    end if  ! Tree / grass IF condition

  end if ! Present and estab IF condition

  !-------------------------
  ! Update LAI and FPC
  if (present(pft)) then

    if (crownarea(pft) > 0.0) then
      lai_ind(pft) = (lm_ind(pft) * sla(pft)) / crownarea(pft)
    else
      lai_ind(pft) = 0.0
    end if

    !------

    fpc_ind(pft)  = (1. - exp(-0.5 * lai_ind(pft)))
    fpc_grid(pft) = crownarea(pft) * nind(pft) * fpc_ind(pft)

  end if

end do

!-------------------------

! if (lprint .and. grid==gprint) write(0,*) lm_ind, rm_ind, sm_ind, hm_ind, nind, &
!                                         present, estab, survive, lai_ind, fpc_ind, fpc_grid

end subroutine establishment

!---------------------------------------------------------------------

end module establishmentmod
