module phenologymod

! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine summerphenology(ndyear,mtemp,dtemp,dphen_t)

! Derive temperature-based phenology for summergreen PFTs

use parametersmod, only : i4,sp,midday
use pftparmod,     only : npft,pftpar,tree,summergreen

integer(i4),                 intent(in)    :: ndyear
real(sp),    dimension(:),   intent(in)    :: mtemp      ! Monthly mean temperature of the current year (12 months) (degC)
real(sp),    dimension(:),   intent(in)    :: dtemp      ! Daily mean temperature (degC)
real(sp),    dimension(:,:), intent(inout) :: dphen_t    ! Temperature based phenology of summergreen (proportion of leaf-on) (fraction)

!-------------------------
! Parameters
real(sp), parameter :: gdd5base = 5.        ! Base temperature for gdd5 calculation
! real(sp), parameter :: ramp     = 1000.     ! Number of GDD5 to attain full leaf cover (from pftpar LPJ)

!-------------------------
! Local variables
integer(i4) :: warmest
integer(i4) :: coldest
integer(i4) :: midsummer
integer(i4) :: firstday

real(sp) :: ramp
real(sp) :: gdd
real(sp) :: leafon
real(sp) :: aphen

integer(i4) :: pft
integer(i4) :: day
integer(i4) :: d

!-------------------------
! Find warmest month, colest month and day of midsummer

warmest = maxloc(mtemp, dim=1)
coldest = minloc(mtemp, dim=1)

midsummer = midday(warmest)

!-------------------------

firstday = midsummer + 1

do while (dtemp(firstday) >= gdd5base .and. firstday /= midsummer)
  firstday = firstday + 1
  if (firstday > ndyear) firstday = 1
end do

!-------------------------

do pft = 1, npft

  if (summergreen(pft)) then

    ramp = pftpar%gdd5ramp(pft)

    if (firstday == midsummer) then

      do day = 1, ndyear
        dphen_t(day,pft) = 1.0
      end do

    else

      gdd = 0.
      leafon = 0.

      day = firstday + 1

      if (day > ndyear) day = 1

      do while (day /= firstday)

        if (dtemp(day) > gdd5base) then

          gdd = gdd + dtemp(day) - gdd5base

          if (ramp > 0.) then
            leafon = min(gdd / ramp, 1.0)
          else
            leafon = 1.
          end if

        end if

        dphen_t(day,pft) = leafon

        day = day + 1

        if (day > ndyear) day = 1

      end do ! Day loop

    end if  ! Firstday == midsummer loop

    !-------------------------
    ! Constrain woody deciduous summergreen phenology to max= 9 months

    if (tree(pft)) then   ! Summergreen trees

      aphen = 0.

      do day = 1, ndyear
        aphen = aphen + dphen_t(day,pft)
      end do

      if (aphen > 210.) then

        do d  = midday(coldest), midday(coldest)+75

          if (d < ndyear) day = d
          if (d > ndyear) day = d - ndyear

          dphen_t(day,pft) = 0.

        end do

        !------

        do d = midday(coldest)-75, midday(coldest)

          if (d > 1) then
            day = d
          else
            day = ndyear + d
          end if

          dphen_t(day,pft) = 0.

        end do

      end if

    end if ! Tree IF condition

  !-------------------------

  else ! Non-summergreen taxa

  !-------------------------

  do day = 1, ndyear
    dphen_t(day,pft) = 1.0
  end do

  end if ! Summergreen IF condition

end do  ! PFT loop


!-------------------------

! if (lprint .and. grid == gprint) write(0,*) dphen_t(:,3)


end subroutine summerphenology




end module phenologymod
