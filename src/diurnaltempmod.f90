module diurnaltempmod

! Copied from arve-dgvm https://github.com/jedokaplan/arve-dgvm_original.git
! Edited by Leo O Lai, HKU, 2021

implicit none

public :: diurnaltemp

!-------
contains

! Formulas used herein are adapted from Cesaraccio, C. Spano, D., Pierpaolo, D., Snyder R. Int. J. Biometerol.
! 2001 (45) 161-169
! NOTE there is an error in the paper,  7a should have a sin after the alpha.
!--

!-----------------------------------------------------------

subroutine diurnaltemp(day,lat,tmin,tmin_n,tmax,tmean,tday,tnight,dayl,  &
                       sunrise,sunset,sunrise_n,dayhour,nighthour)

use parametersmod, only : i4,sp,dp,pi,d2r,daysec
use daylengthmod,  only : daylength

implicit none

integer(i4), intent(in)  :: day
real(dp),    intent(in)  :: lat
real(sp),    intent(in)  :: tmin            ! 24 hour mean minimum temperature (degC)
real(sp),    intent(in)  :: tmin_n          ! Next day minimum temperature
real(sp),    intent(in)  :: tmax            ! 24 hour mean maximum temperature (degC)
real(sp),    intent(out) :: tmean           ! 24 hour mean temperature (degC)
real(sp),    intent(out) :: tday            ! Mean daytime temperature (degC)
real(sp),    intent(out) :: tnight          ! Mean nighttime temperature (degC)
real(sp),    intent(out) :: dayl
integer(i4), intent(out) :: sunrise
integer(i4), intent(out) :: sunset
integer(i4), intent(out) :: sunrise_n
integer(i4), intent(out) :: dayhour
integer(i4), intent(out) :: nighthour

! Parameter
real(dp), parameter :: tfpk = 1./6.                 !delay between solar noon and peak temperature (fraction)

! Local variables
real(dp), dimension(2) :: dayll
real(dp) :: hdl,hdlnext             !half day length (sec) (this day and next)
real(dp) :: ti1                     !midnight till sunup
real(dp) :: ti                      !sundown till midnight
real(dp) :: tam                     !daytime till noon
real(dp) :: tpm                     !daytime post noon till sundown
real(dp) :: morn                    !sunrise - peakt
real(dp) :: sunrise_next            !relative to solar noon (rad)

real(dp) :: sunrse         ! Time of (relative to solar noon (rad))
real(dp) :: sundwn         ! Time of (relative to solar noon (rad))
real(dp) :: peakt          ! Time of peak temperature
real(dp) :: t0             ! Temperature at the sunset hour
real(dp) :: r              ! Difference between max temp and temp at sunset (C)
real(dp) :: a              ! Amplitude of temp change (max - min) in that day (C)
real(dp) :: b              ! Eqn 10 from Cesaraccio paper (ref below)
real(dp) :: hni            ! Length of the night (hrs)

integer(i4) :: dayhour_n

!-------------------------
! Calculate current day and next day length
! Subroutine calc dayl of next day, so -1 will give current day
call daylength(day-1, lat, dayll(1))
call daylength(day, lat, dayll(2))

! For seperating daytime and nighttime data points: (Leo Lai May 2019)
! Define daylength of polar night to be 1 hour
! Define daylength of polar day to be 23 hours
if (dayll(1) >= daysec - 3600.) dayll(1) = daysec - 3600.
if (dayll(2) >= daysec - 3600.) dayll(2) = daysec - 3600.

if (dayll(1) < 3600.) dayll(1) = 3600.
if (dayll(2) < 3600.) dayll(2) = 3600.

!these are prep for other calculations below
!tmin_n = minimum temperature of the following day
t0 = tmax - 0.39 * (tmax - tmin_n)
a = tmax - tmin
r = tmax - t0

!find sunrise and sundown
hdl = 0.5 * dayll(1) / 3600.  !hrs
sunrse = 12. - hdl  !hrs, fixes time from noon till sun up by using half the day length
sundwn = 12. + hdl

!find next days sunrise
hdlnext = 0.5 * dayll(2) / 3600.
sunrise_next = 12 - hdlnext

!this gets the time of peak temperature
peakt = 12. + 2. * hdl * tfpk

! if (dayl(1) <= 3600.) then
!   !FLAG this could be a bit of an ugly way of doing this... it could bias the weather to be mostly
!   !hot since the day (which is the longest) will be set to the max temp of the day. The main reason for
!   !this problem is that for gridded data we do not have the mean daily temp calced out in weathergen
!   !like we do for the max and min. Will think of better ways to do this -JM Oct 29 08
!       !FLAG check on this again after we fix weathergen JM Jan 10 09
!   tday = tmax
!   tnight = tmin

! if (dayl(1) == dayl(2) .AND. dayl(1) /= 0.) print *, dayl

if (dayll(1) < daysec) then !daylength is less than 24 hours

        !has a night and a day, calculate night first

          !find the length of the night
          hni = (43200. - 0.5 * dayll(1) + 43200. - 0.5 * dayll(2))  / 3600. !hrs

          b = (tmin_n - t0) / sqrt(hni) !eqn 10 from Cesaraccio paper

          ti  = t0 * sundwn                                            !sundown
          ti1 = t0 * (sundwn + hni) + 2./3. * b * (hni**(3./2.))   !sunrise (next morn)

          tnight = (ti1 - ti) / hni

          if (dayll(1) > 0.) then  !regular night and day

                    !morning integral (ti is at sunrise, ti1 is at temperature peak time)
                    morn = sunrse - peakt

                    ti  = (tmin * sunrse) + (1. / pi) * (2. * a * morn)
                    ti1 = (tmin * peakt)   + (1. / pi) * (2. * a * morn * cos(pi/2. * (sunrse - peakt) / morn))

                    tam = (ti1 - ti) / (-morn)

                    !afternoon integral (ti is at temperature peak time, ti1 is at sundown)
                    ti  = t0 * peakt   - (1. / pi) * 8. * r * cos(pi / 8. * (-4.))
                    ti1 = t0 * sundwn - (1. / pi) * 8. * r * cos(pi / 8. * (peakt - sundwn - 4.))

                    tpm = (ti1 - ti) / (sundwn - peakt)

                    tday = (tam + tpm) / 2.

          ! else
          !
          !       tday = tnight      !only night, day = night

          end if

end if

! COMMENTED OUT BY Leo Lai (May 2019) because polar day and polar night is redefined above
! else !no night, only day
!
!           !morning integral (ti is at sunrise, ti1 is at temperature peak time)
!           morn = sunrise - peakt
!
!           ti  = tmin * sunrise + 1. / pi * (2. * a * morn)
!           ti1 = tmin * peakt   + 1. / pi * (2. * a * morn * cos(pi / 2. * (sunrise - peakt) / morn))
!
!           tam = (ti1 - ti) / (-morn)
!
!           !afternoon integral (t10 is at temperature peak time, ti1 is at sundown)
!           ti  = t0 * peakt   - 1. / pi * 8. * r * cos(pi / 8. * (-4.))
!           ti1 = t0 * sundown - 1. / pi * 8. * r * cos(pi / 8. * (peakt - sundown - 4.))
!
!
!           tpm = (ti1 - ti) / (sundown - peakt)
!
!           tday = (tam + tpm) / 2.
!
!           tnight = tday
!
! end if

! Calculate mean temperature from tday and tnight to stored into 'dayvars'
tmean = tday * (dayll(1) / daysec) + tnight * ((daysec - dayll(1)) / daysec)

! Claculate and save daylength variables to dayvars
dayl = dayll(1) / 3600.      ! second to hour
! dayl_n = dayl(2) / 3600.    ! second to hour

! Find number of hours with sunlight on current and next day (closest integer for indexing in sub-daily processes)
dayhour   = ceiling(dayll(1) / 3600.)
dayhour_n = ceiling(dayll(2) / 3600.)

! Find current day sunrise and sunset hour
! Variable adapted for indexing on a 24 hourly array such that
!     sunrise = the first hour at which sunlight is observed (inclusive) --> same for sunrise_n
!     sunset = the first hour at which there is NO sunlight (i.e. if sunset = 18:00, meaning light is last observed in 17:00 hour interval)
! NOTE: for indexing, day == (sunrise:sunset-1) and night == (sunset:24) + (1:sunrise_n-1)

if (mod(dayhour,2) /= 0) then

  sunrise = 12 - (dayhour-1) / 2
  sunset  = 12 + (dayhour-1) / 2 + 1

else

  sunrise = 12 - (dayhour/2) + 1
  sunset  = 12 + (dayhour/2) + 1

end if

!---

if (mod(dayhour_n,2) /= 0) then

  sunrise_n = 12 - (dayhour_n-1) / 2

else

  sunrise_n = 12 - (dayhour_n/2) + 1

end if

!---

nighthour = (24 - sunset + 1) + (sunrise_n - 1)


! if (lprint .and. grid==gprint) &
!   print*, dayhour, nighthour, dayvars(grid,day)%sunset-dayvars(grid,day)%sunrise

! if (dayvars(grid,day)%sunset-dayvars(grid,day)%sunrise /= dayhour) &
!   print *, dayhour, nighthour, dayvars(grid,day)%sunrise, dayvars(grid,day)%sunrise_n

end subroutine diurnaltemp

!-----------------------------------------------------------

subroutine humidity(tmean,tdew,rhum)

! Calculate relative humidity from dew point temperature
! From Lawrence (2005) The relationship between relative humidity and the dewpoint temperature. A. MetSoc (https://doi.org/10.1175/BAMS-86-2-225)
! Equation (11)

use parametersmod, only : sp,dp

implicit none

real(sp), intent(in)  :: tmean
real(sp), intent(in)  :: tdew
real(sp), intent(out) :: rhum

real(sp) :: tmean_K
real(sp) :: tdew_K

real(dp), parameter :: L  = 2257000
real(dp), parameter :: Rw = 461.5

!------

tmean_K = tmean + 273.15
tdew_K  = tdew + 273.15

!------

rhum = 100. * exp(((1. - (tmean_K / tdew_K)) * (L / Rw)) / tmean_K)

if (rhum > 100.) rhum = 100.    ! When temp is lower than dew point


end subroutine humidity


end module diurnaltempmod
