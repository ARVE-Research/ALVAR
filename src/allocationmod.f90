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

use pftparmod,    only : npft,pftpar,tree
use statevarsmod, only : ndyear,dayvars,soilvars,gppvars,vegvars,topovars,lprint,gprint,gridlon,gridlat

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
real(sp), parameter :: nind_min         = 1.e-10      !minimum individual density for persistence of PFT (indiv/m2)
real(sp), parameter :: eps              = 1.e-6       !epsilon parameter (?)
integer, parameter :: nseg = 20
real(sp), parameter :: xacc =  0.1     !x-axis precision threshold for the allocation solution
real(dp), parameter :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution

real(sp), parameter :: crownarea_max = 15.0

!-------------------------
! Pointer variables
logical,  pointer, dimension(:) :: present
logical,  pointer, dimension(:) :: estab
logical,  pointer, dimension(:) :: survive
real(sp), pointer, dimension(:) :: dwscal
real(sp), pointer, dimension(:) :: fpc_grid
real(sp), pointer, dimension(:) :: fpc_ind
real(sp), pointer, dimension(:) :: fpc_inc
real(sp), pointer, dimension(:) :: nind
real(sp), pointer, dimension(:) :: sla
real(sp), pointer, dimension(:) :: stemdiam
real(sp), pointer, dimension(:) :: height
real(sp), pointer, dimension(:) :: crownarea
real(sp), pointer, dimension(:) :: lai_ind
real(sp), pointer, dimension(:) :: hm_ind
real(sp), pointer, dimension(:) :: lm_ind
real(sp), pointer, dimension(:) :: sm_ind
real(sp), pointer, dimension(:) :: rm_ind
real(sp), pointer, dimension(:) :: npp

real(sp), pointer, dimension(:) :: litter_ag_fast    ! Fast above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_ag_slow    ! Slow above ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: litter_bg         ! Below ground litter pool (gC m-2)
real(sp), pointer, dimension(:) :: bm_inc

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
real(sp) :: ir_max

real(sp) :: lminc_ind
real(sp) :: rminc_ind
real(sp) :: sminc_ind

real(sp) :: fpc_grid_old

real(sp), dimension(4) :: treecarbon

integer :: i
integer :: pft


!-------------------------

npp       => gppvars(grid,day)%npp
dwscal    => gppvars(grid,day)%dwscal

present   => vegvars(grid)%present
estab     => vegvars(grid)%estab
survive   => vegvars(grid)%survive
fpc_grid  => vegvars(grid)%fpc_grid
fpc_ind   => vegvars(grid)%fpc_ind
fpc_inc   => vegvars(grid)%fpc_inc
nind      => vegvars(grid)%nind
sla       => vegvars(grid)%sla
stemdiam  => vegvars(grid)%stemdiam
height    => vegvars(grid)%height
crownarea => vegvars(grid)%crownarea
lai_ind   => vegvars(grid)%lai_ind
hm_ind    => vegvars(grid)%hm_ind
lm_ind    => vegvars(grid)%lm_ind
sm_ind    => vegvars(grid)%sm_ind
rm_ind    => vegvars(grid)%rm_ind

litter_ag_fast => vegvars(grid)%litter_ag_fast
litter_ag_slow => vegvars(grid)%litter_ag_slow
litter_bg      => vegvars(grid)%litter_bg
bm_inc         => vegvars(grid)%bm_inc

!-------------------------

do pft = 1, npft

  if (present(pft)) then

    if (npp(pft) == -9999.) return

    lm = lm_ind(pft)
    sm = sm_ind(pft)
    hm = hm_ind(pft)
    rm = rm_ind(pft)

    bm_inc(pft) = sum(gppvars(grid,:)%npp(pft))

    bm_inc_ind = bm_inc(pft) / nind(pft)

    bm_inc_ind = 0.9 * bm_inc_ind         ! 10% reproduction cost inserted here for now (Leo Lai Oct 2021)

    ! if (bm_inc_ind < 1.e-6) write(0,*) 'pft: ', pft, 'needs to be killed'

    ! if(lprint .and. grid == gprint) print *, sum(vegvars(grid,:)%gpp), sum(vegvars(grid,:)%npp0), sum(vegvars(grid,:)%aresp)

    ! awscal = sum(vegvars(grid,:)%dwscal) / ndyear

    ir_max = pftpar(16,pft)

    awscal = dwscal(pft)

    lm2rm = max(ir_max * awscal, 0.1)

    if (tree(pft)) then

      !====================
      ! TREE ALLOCATION
      !====================

      lm1 = latosa * sm / (wooddens * height(pft) * sla(pft))     ! Allometric leaf mass requirement

      lminc_ind_min = lm1 - lm        ! Minimum leaf mass increment based on current allometry

      ! if (lprint .and. grid==gprint) print *, lm1, lm

      ! Calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
      ! May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.

      rminc_ind_min = lm1 / lm2rm - rm

      !-------------------------

      if (rminc_ind_min > 0. .and. lminc_ind_min > 0. .and. rminc_ind_min + lminc_ind_min <= bm_inc_ind) then

        if (lprint .and. grid==gprint) write(0,*) 'Normal allocation'

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

          fx1 = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla(pft),x1)

          !Find approximate location of leftmost root on the interval (x1,x2).
          !Subdivide (x1,x2) into nseg equal segments seeking change in sign of f(xmid) relative to f(x1).

          fmid = fx1
          xmid = x1

          i = 1

          do

            xmid = xmid + dx

            fmid = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla(pft),xmid)

            if (fmid * fx1 <= 0. .or. xmid >= x2) exit  !sign has changed or we are over the upper bound

            if (i > 20) write(0,*)'first alloc loop flag',i,fmid*fx1,xmid,x1,x2,dx,bm_inc_ind
            if (i > 50) stop 'Too many iterations allocmod'

            i = i + 1

          end do

          !the interval that brackets zero in f(x) becomes the new bounds for the root search

          x1 = xmid - dx
          x2 = xmid

          !Apply bisection method to find root on the new interval (x1,x2)

          fx1 = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla(pft),x1)

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

            fmid = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla(pft),xmid)

            if (fmid * sign <= 0.) rtbis = xmid

            if (dx < xacc .or. abs(fmid) <= yacc) exit

            if (i > 20) write(0,*)'second alloc loop flag',i,dx,abs(fmid),'pft:',pft,fpc_grid,nind,gridlon(grid),gridlat(grid)
            if (i > 50) stop 'Too many iterations allocmod'

            i = i + 1

          end do

          !Now rtbis contains numerical solution for lminc_ind given eqn (22)

          lminc_ind = rtbis

        end if  !x2-x1 block

        !Calculate increments in other compartments using allometry relationships

        rminc_ind = (lm + lminc_ind) / lm2rm - rm       !eqn (9)

        sminc_ind = bm_inc_ind - lminc_ind - rminc_ind  !eqn (1)

      !-------------------------

      else

      !-------------------------

        ! if (lprint .and. grid==gprint) write(0,*) 'Abnormal allocation'

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

        sminc_ind = (lm + lminc_ind) * sla(pft) / latosa * wooddens * height(pft) - sm  !eqn (35)

        ! if (lprint .and. grid==gprint) print *, sminc_ind

        !Convert killed sapwood to heartwood

        hm = hm + abs(sminc_ind)

        !write(stdout,*)'abnormal case',i,lminc_ind,rminc_ind,lminc_ind+rminc_ind,bm_inc_ind,sminc_ind

        ! if (lprint .and. grid==gprint .and. year > 1) then
        !   write(0,*) 'Abnormal allocation'
        !   print *, year, day0, lminc_ind, rminc_ind, sminc_ind, bm_inc_ind, &
        !           lm_ind, rm_ind, sm_ind, hm_ind, nind, &
        !           height, stemdiam, crownarea, fpc_grid
        !   ! stop
        ! end if


      end if  ! Normal/abnormal allocation

      !-------------------------
      ! Increment C compartments

      lm_ind(pft) = lm + lminc_ind
      rm_ind(pft) = rm + rminc_ind
      sm_ind(pft) = sm + sminc_ind
      hm_ind(pft) = hm

      !-------------------------
      ! Calculate new height, diameter and crown area

      if (lm_ind(pft) > 0.) then

        sap_xsa = lm_ind(pft) * sla(pft) / latosa  !eqn (5)

        height(pft) = sm_ind(pft) / sap_xsa / wooddens

        stemdiam(pft)    = (height(pft) / allom2) ** (1. / allom3)                  !eqn (C)

        crownarea(pft) = min(allom1 * (stemdiam(pft)**reinickerp), crownarea_max)  !eqn (D)

      end if

    else ! Grass allocation IF condition

      !====================
      ! GRASS ALLOCATION
      !====================

      ! Distribute this year's production among leaves and fine roots according to leaf to rootmass ratio [eqn (33)] (see below)
      ! Relocation of C from one compartment to the other not allowed: negative increment in either compartment transferred to litter
      ! but the total negative amount cannot be more than the existing pool plus the increment

      lminc_ind = (bm_inc_ind - lm / lm2rm + rm) / (1. + 1. / lm2rm)

      rminc_ind = bm_inc_ind - lminc_ind

      if (lminc_ind > 0.) then

        if (rminc_ind < 0.) then  ! Negative allocation to grass root mass

          if (rminc_ind + rm < 0.) rminc_ind = -rm  ! Cannot be more than the grass root mass that is actually present

          ! Add killed grass roots to below-ground litter

          litter_bg(pft) = litter_bg(pft) + abs(rminc_ind) * nind(pft)

        end if

      else

        ! Negative allocation to grass leaf mass

        rminc_ind = bm_inc_ind
        lminc_ind = lm2rm * (rm + rminc_ind) - lm

        if (lminc_ind > 0.) lminc_ind = -lm  !cannot be more than the grass leaf mass that is actually present

        ! Add killed grass leaf mass to litter

        litter_ag_fast(pft) = litter_ag_fast(pft) + abs(lminc_ind) * nind(pft)

      end if

      ! Increment grass C compartments

      lm_ind(pft) = lm + lminc_ind
      rm_ind(pft) = rm + rminc_ind


    end if ! Tree / grass IF condition

    !-------------------------
    ! Update LAI and FPC

    if (crownarea(pft) > 0.) then
      lai_ind(pft) = (lm_ind(pft) * sla(pft)) / crownarea(pft)
    else
      lai_ind(pft) = 0.
    end if

    fpc_grid_old  = fpc_grid(pft)
    fpc_ind(pft)  = 1. - exp(-0.5 * lai_ind(pft))
    fpc_grid(pft) = crownarea(pft) * nind(pft) * fpc_ind(pft)
    fpc_inc(pft)  = max(fpc_grid(pft) - fpc_grid_old, 0.)

  end if ! Present / absent IF condition

end do ! PFT loop


!-------------------------
! Check validity of allocation and correct (from lpjmod.f90 in LPJ-LMFire)
! Heartwood can be zero, but all other pools have to be positive to have valid allometry

do pft = 1, npft

  if (.not.tree(pft)) cycle

  treecarbon(1) = lm_ind(pft)
  treecarbon(2) = rm_ind(pft)
  treecarbon(3) = sm_ind(pft)
  treecarbon(4) = hm_ind(pft)

  if (any(treecarbon(1:3) <= 0.) .and. (sum(treecarbon) > 0. .or. nind(pft) > 0.)) then

    ! write(0,*) 'Invalid allometry after allocation in year ', year, 'for grid ', grid, 'NPP:', bm_inc_ind

    litter_ag_fast(pft) = litter_ag_fast(pft) + nind(pft) * lm_ind(pft)
    litter_ag_slow(pft) = litter_ag_slow(pft) + nind(pft) * (sm_ind(pft) + hm_ind(pft))
    litter_bg(pft)      = litter_bg(pft)      + nind(pft) * rm_ind(pft)

    present(pft)  = .false.

    nind(pft)   = 0.
    lm_ind(pft) = 0.
    rm_ind(pft) = 0.
    sm_ind(pft) = 0.
    hm_ind(pft) = 0.

    ! fpc_grid(pft) = 0.

  end if

end do

!-------------------------


if (lprint .and. grid==gprint) write(0,*) year, day, lminc_ind, rminc_ind, sminc_ind, bm_inc_ind, &
                                        lm_ind(4), rm_ind(4), sm_ind(4), hm_ind(4), nind(4), &
                                        height(4), stemdiam(4), crownarea(4), lai_ind(4), &
                                        sum(dwscal)/npft, sum(dayvars(grid,:)%prec)

if (lprint .and. grid==gprint) write(0,*) "bm_inc: ", bm_inc, "bm_inc_ind: ", 0.9*  bm_inc / nind, 'fpc_grid:', fpc_grid

! if (any(vegvars(grid)%fpc_grid > 1.0)) write(0,*) year, vegvars(grid)%fpc_grid
!
! if (lprint .and. grid==gprint) print *, sum(vegvars(grid,:)%gpp), &
!                                         sum(vegvars(grid,:)%aresp), &
!                                         sum(vegvars(grid,:)%npp)


! if (height > 20.) print *, gridlon(grid), gridlat(grid)









end subroutine allocation

!---------------------------------------------------------

real(sp) function root(lm,sm,hm,rm,inc,lm2rm,sla,x)

implicit none

!parameters
real(sp), parameter :: pi    = 3.1415926535
real(sp), parameter :: pi4    = pi / 4.
real(sp), parameter :: allom1 = 100.
real(sp), parameter :: allom2 = 40.
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
