module allocationmod

! Module to calculate vegation carbon allocation
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine allocation(year,grid,day)

use statevarsmod, only : ndyear,dayvars,soilvars,vegvars,topovars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: pi    = 3.1415926535
real(sp), parameter :: allom1 = 100
real(sp), parameter :: allom2 = 40
real(sp), parameter :: allom3 = 0.5
real(sp), parameter :: allom4 = 0.3
real(sp), parameter :: reinickerp = 1.6
real(sp), parameter :: latosa     = 8.e3
real(sp), parameter :: wooddens   = 2.e5
real(sp), parameter :: aprec_min_estab  = 100.        !minimum annual precipitation for establishment (mm)
real(sp), parameter :: estab_max        = 0.24        !maximum sapling establishment rate (indiv/m2)
real(sp), parameter :: nind_min         = 1.e-10      !minimum individual density for persistence of PFT (indiv/m2)
real(sp), parameter :: eps              = 1.e-6       !epsilon parameter (?)
real(sp), parameter :: ir_max              = 1.0      ! Leaf-to-root ratio under nonwater stressed conditions
integer, parameter :: nseg = 20
real(sp), parameter :: xacc =  0.1     !x-axis precision threshold for the allocation solution
real(dp), parameter :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution

real(sp), parameter :: crownarea_max = 15.0

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
real(sp), pointer :: npp
real(sp), pointer :: gpp

!-------------------------
! Local variables
real(sp) :: lm
real(sp) :: rm
real(sp) :: sm
real(sp) :: hm
real(sp) :: bm_inc_ind
real(sp) :: awscal
real(sp) :: lm2rm
real(sp) :: lminc_ind_min
real(sp) :: lm1
real(sp) :: rminc_ind_min
logical :: normal
real(sp) :: x1
real(sp) :: x2
real(sp) :: dx
real(sp) :: fx1
real(sp) :: fmid
real(sp) :: xmid
real(sp) :: sign
real(sp) :: rtbis
real(sp) :: sap_xsa

real(sp) :: lminc_ind
real(sp) :: rminc_ind
real(sp) :: sminc_ind

real(sp) :: fpc_grid_old


integer :: i


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
npp => vegvars(grid,day)%npp0
gpp => vegvars(grid,day)%gpp0

!-------------------------

if (present) then

  if (npp == -9999.) vegvars(grid,:)%npp0 = 0.0

  lm = lm_ind
  sm = sm_ind
  hm = hm_ind
  rm = rm_ind

  bm_inc_ind = sum(vegvars(grid,:)%npp0) / nind

  bm_inc_ind = 0.9 * bm_inc_ind         ! 10% reproduction cost inserted here for now (Leo Lai Oct 2021)

  ! if(lprint .and. grid == gprint) print *, sum(vegvars(grid,:)%gpp), sum(vegvars(grid,:)%npp0), sum(vegvars(grid,:)%aresp)

  awscal = sum(vegvars(grid,:)%dwscal) / ndyear

  lm2rm = max(ir_max * awscal, 0.1)

  !------
  ! Tree allocation

  lm1 = latosa * sm / (wooddens * height * sla)     ! Allometric leaf mass requirement

  lminc_ind_min = lm1 - lm        ! Minimum leaf mass increment based on current allometry

  ! print *, lm1, lm

  ! Calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
  ! May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.

  rminc_ind_min = lm1 / lm2rm - rm


  if (rminc_ind_min > 0. .and. lminc_ind_min > 0. .and. rminc_ind_min + lminc_ind_min <= bm_inc_ind) then

    if (lprint .and. grid==gprint) print *, 'Normal allocation'

    !Normal allocation (positive increment to all living C compartments)

    normal = .true.

    !Calculation of leaf mass increment (lminc_ind) that satisfies Eqn (22)
    !Since this is normal allocation, we set the lower bound for the leafmass allocation (x1)
    !to its allometric minimum, because it should be able to be fulfilled, i.e.:

    x1 = lminc_ind_min
    x2 = (bm_inc_ind - (lm / lm2rm - rm)) / (1. + 1. / lm2rm)

    dx = x2 - x1

    if (dx < 0.01) then

      !there seems to be rare cases where lminc_ind_min (x1) is almost equal to x2. In this case,
      !assume that the leafmass increment is equal to the midpoint between the values and skip
      !the root finding procedure

      lminc_ind = x1 + 0.5 * dx

    else

      !Find a root for non-negative lminc_ind, rminc_ind and sminc_ind using Bisection Method (Press et al 1986, p 346)
      !There should be exactly one solution (no proof presented, but Steve has managed one).

      dx = dx / real(nseg)

      !evaluate f(x1) = LHS of eqn (22) at x1

      fx1 = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla,x1)

      !Find approximate location of leftmost root on the interval (x1,x2).
      !Subdivide (x1,x2) into nseg equal segments seeking change in sign of f(xmid) relative to f(x1).

      fmid = fx1
      xmid = x1

      i = 1

      do

        xmid = xmid + dx

        fmid = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla,xmid)

        if (fmid * fx1 <= 0. .or. xmid >= x2) exit  !sign has changed or we are over the upper bound

        if (i > 20) write(0,*)'first alloc loop flag',i,fmid*fx1,xmid,x1,x2,dx,bm_inc_ind
        if (i > 50) stop 'Too many iterations allocmod'

        i = i + 1

      end do

      !the interval that brackets zero in f(x) becomes the new bounds for the root search

      x1 = xmid - dx
      x2 = xmid

      !Apply bisection method to find root on the new interval (x1,x2)

      fx1 = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla,x1)

      if (fx1 >= 0.) then
        sign = -1.
      else
        sign =  1.
      end if

      rtbis = x1
      dx    = x2 - x1

      !Bisection loop: search iterates on value of xmid until xmid lies within xacc of the root,
      !i.e. until |xmid-x| < xacc where f(x) = 0. the final value of xmid with be the leafmass increment

      i = 1

      do

        dx   = 0.5 * dx
        xmid = rtbis + dx

        !calculate fmid = f(xmid) [eqn (22)]

        fmid = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla,xmid)

        if (fmid * sign <= 0.) rtbis = xmid

        if (dx < xacc .or. abs(fmid) <= yacc) exit

        if (i > 20) write(0,*)'second alloc loop flag',i,dx,abs(fmid)
        if (i > 50) stop 'Too many iterations allocmod'

        i = i + 1

      end do

      !Now rtbis contains numerical solution for lminc_ind given eqn (22)

      lminc_ind = rtbis

    end if  !x2-x1 block

    !Calculate increments in other compartments using allometry relationships

    rminc_ind = (lm + lminc_ind) / lm2rm - rm       !eqn (9)

    sminc_ind = bm_inc_ind - lminc_ind - rminc_ind  !eqn (1)

  else

    if (lprint .and. grid==gprint) print *, 'Abnormal allocation'

    !Abnormal allocation: reduction in some C compartment(s) to satisfy allometry

    normal = .false.

    !Attempt to distribute this year's production among leaves and roots only

    lminc_ind = (bm_inc_ind - lm / lm2rm + rm) / (1. + 1. / lm2rm)  !eqn (33)

    ! if (lprint .and. grid==gprint) print *, lminc_ind

    if (lminc_ind > 0.) then

      !Positive allocation to leafmass

      rminc_ind = bm_inc_ind - lminc_ind  !eqn (31)

      ! if (lprint .and. grid==gprint) print *, rminc_ind

      !Add killed roots (if any) to below-ground litter

      if (rminc_ind < 0.) then

        lminc_ind = bm_inc_ind
        rminc_ind = (lm + lminc_ind) / lm2rm - rm

        ! litter_bg(pft,1) = litter_bg(pft,1) + abs(rminc_ind) * nind

      end if

      i = 1

    else

      !Negative allocation to leaf mass

      rminc_ind = bm_inc_ind
      lminc_ind = (rm + rminc_ind) * lm2rm - lm  !from eqn (9)

      !Add killed leaves to litter

      ! litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + abs(lminc_ind) * nind

      i = 2

    end if

    !Calculate sminc_ind (must be negative)

    sminc_ind = (lm + lminc_ind) * sla / latosa * wooddens * height - sm  !eqn (35)

    ! if (lprint .and. grid==gprint) print *, sminc_ind

    !Convert killed sapwood to heartwood

    hm = hm + abs(sminc_ind)

    !write(stdout,*)'abnormal case',i,lminc_ind,rminc_ind,lminc_ind+rminc_ind,bm_inc_ind,sminc_ind

  end if  !normal/abnormal allocation

  !################################
  ! Testing site -> simple assumption of allocation ratio (arbitarily derived based on Sampling mass- Leo)

  ! lminc_ind = 0.25 * bm_inc_ind
  ! rminc_ind = lminc_ind
  ! sminc_ind = bm_inc_ind - lminc_ind - rminc_ind

  !################################


  !Increment C compartments

  lm_ind = lm + lminc_ind
  rm_ind = rm + rminc_ind
  sm_ind = sm + sminc_ind
  hm_ind = hm

  !Calculate new height, diameter and crown area

  if (lm_ind> 0.) then

    sap_xsa = lm_ind * sla / latosa  !eqn (5)

    height = sm_ind / sap_xsa / wooddens

    stemdiam    = (height / allom2) ** (1. / allom3)                  !eqn (C)

    crownarea = min(allom1 * (stemdiam**reinickerp), crownarea_max)  !eqn (D)

  end if

  !Update LAI and FPC

  if (crownarea > 0.) then
    lai_ind = (lm_ind * sla) / crownarea
  else
    lai_ind = 0.
  end if

  fpc_grid_old  = fpc_grid
  fpc_ind       = 1. - exp(-0.5 * lai_ind)
  fpc_grid = crownarea * nind * fpc_ind
  fpc_inc  = max(fpc_grid - fpc_grid_old, 0.)

end if



if (lprint .and. grid==gprint) print *, lminc_ind, rminc_ind, sminc_ind, bm_inc_ind, &
                                        lm_ind, rm_ind, sm_ind, hm_ind, nind, &
                                        height, stemdiam, crownarea, lai_ind, fpc_grid











end subroutine allocation

!---------------------------------------------------------

real(sp) function root(lm,sm,hm,rm,inc,lm2rm,sla,x)

implicit none

!parameters
real(sp), parameter :: pi    = 3.1415926535
real(sp), parameter :: pi4    = pi/4
real(sp), parameter :: allom1 = 100
real(sp), parameter :: allom2 = 40
real(sp), parameter :: allom3 = 0.5
real(sp), parameter :: allom4 = 0.3
real(sp), parameter :: reinickerp = 1.6
real(sp), parameter :: latosa     = 8.e3
real(sp), parameter :: wooddens   = 2.e5
real(sp), parameter :: a1  = 2. / allom3
real(sp), parameter :: a2  = 1. + a1
real(sp), parameter :: a3  = allom2**a1

!arguments

real(sp), intent(in) :: sla    !specific leaf area
real(sp), intent(in) :: lm2rm  !leaf mass to root mass ratio
real(sp), intent(in) :: lm     !individual leaf mass
real(sp), intent(in) :: sm     !individual sapwood mass
real(sp), intent(in) :: hm     !individual heartwood mass
real(sp), intent(in) :: rm     !individual root mass
real(sp), intent(in) :: inc    !individual biomass increment
real(sp), intent(in) :: x      !leafmass allocation amount as input

!---

root = a3 * ((sm + inc - x - ((lm + x) / lm2rm) + rm + hm) / wooddens) / pi4 - &
            ((sm + inc - x - ((lm + x) / lm2rm) + rm) / ((lm + x) * sla * wooddens / latosa))**a2

end function root

!---------------------------------------------------------





end module allocationmod
