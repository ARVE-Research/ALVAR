module fireindexmod

use parametersmod, only : i2,i4,sp,dp

implicit none


! Global variable
integer(i4), allocatable, dimension(:) :: dryd       ! Number of dry days preceding current day / day since rain
real(sp), allocatable, dimension(:) :: KBDI_0        ! Previous day KBDI
real(sp), allocatable, dimension(:) :: prec_acc      ! Rainfall accumulation of the current wet period (mm)


contains

!---------------------------------------------------------------------

subroutine fireindex(rank,year,grid,day)

! Subroutine to calculate McArthur's fire-danger meters
! Equations from Noble et al. (1980) McArthur's fire danger meters expressed as equations https://doi.org/10.1111/j.1442-9993.1980.tb01243.x
! Coded by Leo Lai as part of ALVAR model extension (Jul 2021)

use metvarsmod, only : dayvars,genvars,ndyear,cnt,topovars

integer(i4), intent(in) :: rank
integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

! Pointers
real(sp), pointer :: DF             ! Drought factor
real(sp), pointer :: KBDI           ! Keetch-Byram Drought Index (mm equivalent)
real(sp), pointer :: ForFireMk5     ! Forest fire danger index Mark 5 meter

! Local variables
real(sp) :: tmean         ! 24 hour mean temperature (degC)
real(sp) :: tmax          ! Daily maximum temperature (degC)
real(sp) :: rhum          ! Relative humidity (%)
real(sp) :: wind          ! 10m windspeed (m s-1)
real(sp) :: prec          ! Daily precipitation (mm)
real(sp) :: aprec         ! Annual precipitation (mm)

real(sp) :: ForFireMk5_2

real(sp) :: sloperad

!----------
! Initialize global variables at first time step and first grid

if (year == 1 .and. grid == 1 .and. day == 1) then

  allocate(KBDI_0(cnt(1)))
  allocate(dryd(cnt(1)))
  allocate(prec_acc(cnt(1)))

end if

!----------

DF         => dayvars(grid,day)%DF
KBDI       => dayvars(grid,day)%KBDI
ForFireMk5 => dayvars(grid,day)%ForFireMk5


sloperad = topovars(grid)%sloperad

tmean = dayvars(grid,day)%tmean * sloperad
tmax  = dayvars(grid,day)%tmax
rhum  = dayvars(grid,day)%rhum
wind  = dayvars(grid,day)%wind
prec  = dayvars(grid,day)%prec

aprec = sum(genvars%pre(5:16))

!----------

if (year == 1 .and. day == 1) then

  KBDI_0(grid) = 20.
  dryd(grid) = 0.
  prec_acc(grid) = 0.

end if

!----------
! Calculation of the Keetch_Bryam Drought Index (KBDI)
! Equations taken from Janis et al. (2002 Intl J. of Wildland Fire, 2002, 11, 281-289 http://climate.geog.udel.edu/~climate/publication_html/Pdf/JJF_IJWF_02.pdf

if (prec > 0) then      ! if wet day

  dryd(grid) = 0                           ! reset consecutive dry day
  prec_acc(grid) = prec_acc(grid) + prec   ! accumulate rainfall for wet day (mm)

else

  dryd(grid) = dryd(grid) + 1              ! add one to consecutive dry day count
  prec_acc(grid) = 0.                      ! reset precipitation accumulation (mm)

end if

!----------

DF = (800. - KBDI_0(grid)) * (0.968 * exp(0.0875 * tmax + 1.5552) - 8.30) * 1.e-3  &
      / (1. + 10.88 * exp(-1.74e-3 * aprec))

!----------

if (prec == 0. .and. tmax <= 6.78) then

  KBDI = KBDI_0(grid)

else if  (prec == 0. .and. tmax > 6.78) then

  KBDI = KBDI_0(grid) + DF

else if (prec > 0 .and. prec_acc(grid) < 5.1) then

  KBDI = KBDI_0(grid) + DF

else if (prec > 0 .and. prec_acc(grid) >= 5.1) then

  DF = (800. - KBDI_0(grid) + 3.937 * prec_acc(grid)) * (0.968 * exp(0.0875 * tmax + 1.5552) - 8.30) * 1.e-3  &
        / (1. + 10.88 * exp(-1.74e-3 * aprec))

  KBDI = KBDI_0(grid) - 3.937 * prec_acc(grid) + DF

  prec_acc = 0.                  ! Reset rainfall accumulation after subtracting once >> Eq. 2d text explanation

end if

! if (rank == 1 .and. grid == 500) print *, KBDI, DF, prec, tmax, dryd(grid), prec_acc(grid), aprec

!----------

KBDI = max(0., KBDI)

KBDI_0(grid) = KBDI

!----------
! Calculation of Fire Danger Index
! Equations from Noble et al. (1980) McArthur's fire danger meters expressed as equations https://doi.org/10.1111/j.1442-9993.1980.tb01243.x

! wind = wind / 10.         ! TEMPORARY fix to wrong scaling factor of input data....

! print *, wind

wind = wind * 3.6         ! Convert to km h-1

if (wind > 200) wind = wind / 3.
! if (wind > 400) wind = wind / 5.

! wind = (real(day) / 365.) * (14 * 3.6)     ! Arbitrary windspeed (km hr-1)

DF = 0.191 * (KBDI + 104.) * (dryd(grid) + 1) ** 1.5 / (3.52 * (dryd(grid) + 1) ** 1.5 + prec - 1.)

ForFireMk5 = 2.0 * exp(-4.50 + 0.987 * log(DF) - 0.0345 * rhum + 0.0338 * tmean + 0.0234 * wind)

ForFireMk5_2 = 1.25 * DF * exp((tmean - rhum) / 30.0 + 0.0234 * wind)       ! Simplified euqation from Noble et al.

if (ForFireMk5 > 120.) ForFireMk5 = 120.

! if (rank == 0 .and. grid == 5000) print *, KBDI, DF, ForFireMk5, ForFireMk5_2, dryd(grid), prec, rhum, tmean, wind

! if (wind > 100.) print *, KBDI, DF, ForFireMk5, ForFireMk5_2, dryd(grid), prec, rhum, tmean, wind
!
! if (ForFireMk5 < 0.) print *, 'WTFFFFFFFFFFFFFFF'


end subroutine fireindex

!---------------------------------------------------------------------

end module fireindexmod
