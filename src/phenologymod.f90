module phenologymod

! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine summerphenology(year,grid,day0)

! Derive temperature-based phenology for summergreen PFTs

use parametersmod, only : midday
use statevarsmod, only : ndyear,genvars,dayvars,vegvars,lprint,gprint

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day0

!-------------------------
! Parameters
real(sp), parameter :: gdd5base = 5.        ! Base temperature for gdd5 calculation
real(sp), parameter :: ramp     = 1000.     ! Number of GDD5 to attain full leaf cover (from pftpar LPJ)

!-------------------------
! Pointer variables
real(sp), pointer, dimension(:) :: dtemp    ! Daily mean temperature (degC)
real(sp), pointer, dimension(:) :: dphen_t  ! Temperature based phenology of summergreen (proportion of leaf-on) (fraction)

!-------------------------
! Local variables
integer(i4) :: warmest
integer(i4) :: coldest
integer(i4) :: midsummer
integer(i4) :: firstday

real(sp) :: gdd
real(sp) :: leafon
real(sp) :: aphen
integer  :: day
integer  :: d

!-------------------------

dtemp   => dayvars(grid,1:ndyear)%tmean
dphen_t => vegvars(grid)%dphen_t(1:ndyear)

!-------------------------
! Find warmest month, colest month and day of midsummer

warmest = maxloc(genvars%tmp(5:16), dim=1)
coldest = minloc(genvars%tmp(5:16), dim=1)

midsummer = midday(warmest)

!-------------------------

firstday = midsummer + 1

do while (dtemp(firstday) >= gdd5base .and. firstday /= midsummer)
  firstday = firstday + 1
  if (firstday > ndyear) firstday = 1
end do

!-------------------------

if (firstday == midsummer) then ! No lead abscission

  dphen_t = 1.0

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

    dphen_t(day) = leafon

    day = day + 1

    if (day > ndyear) day = 1

  end do

end if

!-------------------------
! Constrain woody deciduous (tree) phenology to maximum of 9 months

! aphen = 0.
!
! do day = 1, ndyear
!   aphen = aphen + dphen_t(day)
! end do
!
! if (aphen > 210.) then
!
!   do d = midday(coldest), midday(coldest)+75
!
!     if (d < ndyear) day = d
!     if (d > ndyear) day = d - ndyear
!
!     dphen_t(day) = 0.
!
!   end do
!
!   !------
!
!   do d = midday(coldest)-75, midday(coldest)
!
!     if (d > 1) then
!       day = d
!     else
!       day = ndyear + d
!     end if
!
!     dphen_t(day) = 0.
!
!   end do
!
! end if


!-------------------------

! if (lprint .and. grid == gprint) print *, dphen_t


end subroutine summerphenology




end module phenologymod
