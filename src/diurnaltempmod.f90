module diurnaltempmod

! Copied from arve-dgvm https://github.com/jedokaplan/arve-dgvm_original.git
! Edited by Leo O Lai, HKU, 2021

use parametersmod, only : i2,i4,sp,dp,pi,d2r,daysec
use daylengthmod,  only : daylength

implicit none

public :: diurnaltemp

real(dp) :: sunrise     !time of (relative to solar noon (rad))
real(dp) :: sundown        !time of (relative to solar noon (rad))
real(dp) :: peakt       !time of peak temperature
real(dp) :: t0                !temperature at the sunset hour
real(dp) :: r                !difference between max temp and temp at sunset (C)
real(dp) :: a           !amplitude of temp change (max - min) in that day (C)
real(dp) :: b                !eqn 10 from Cesaraccio paper (ref below)
real(dp) :: hni         !length of the night (hrs)
real(dp) :: tfollow     !the next days minimum temperature (assigned based on grid or point mode)

!-------
contains

! Formulas used herein are adapted from Cesaraccio, C. Spano, D., Pierpaolo, D., Snyder R. Int. J. Biometerol.
! 2001 (45) 161-169
!NOTE there is an error in the paper,  7a should have a sin after the alpha.
!--

!-----------------------------------------------------------

subroutine diurnaltemp(grid,day)

! use statevars,   only : sv,gridded,dayl
use metvarsmod,  only : dayvars,gridlat

!arguments
integer(i4), intent(in) :: grid
integer, intent(in) :: day       ! Julian day (1 to 366)

!pointers
real(sp), pointer :: tmin                       !min temperature (C)
real(sp), pointer :: tmax                       !max temperature (C)
real(sp), pointer :: tnextmin                   !minimum temperature of the following day
real(sp), pointer :: tday                       !daytime temperature (C)
real(sp), pointer :: tnight                     !nighttime temperature (C)
real(sp), pointer :: tmean                      !mean temperature (C)

real(dp) :: lat

!parameter
real(dp), parameter :: tfpk = 1./6.                 !delay between solar noon and peak temperature (fraction)

!local variables
real(dp), dimension(2) :: dayl
real(dp) :: hdl,hdlnext             !half day length (sec) (this day and next)
real(dp) :: ti1                     !midnight till sunup
real(dp) :: ti                      !sundown till midnight
real(dp) :: tam                     !daytime till noon
real(dp) :: tpm                     !daytime post noon till sundown
real(dp) :: morn                    !sunrise - peakt
real(dp) :: sunrise_next            !relative to solar noon (rad)

!point pointers to global data
tmin     => dayvars(grid,day)%tmin
tmax     => dayvars(grid,day)%tmax
tmean    => dayvars(grid,day)%tmean
tnextmin => dayvars(grid,day+1)%tmin
tday     => dayvars(grid,day)%tday
tnight   => dayvars(grid,day)%tnight

!------
!initial assignments

tfollow = tnextmin

lat = gridlat(grid)

! Calculate current day and next day length
! Subroutine calc dayl of next day, so -1 will give current day
call daylength(day-1, lat, dayl(1))
call daylength(day, lat, dayl(2))

!these are prep for other calculations below
!tfollow = minimum temperature of the following day
t0 = tmax - 0.39 * (tmax - tfollow)
a = tmax - tmin
r = tmax - t0

!find sunrise and sundown
hdl = 0.5 * dayl(1) / 3600.  !hrs
sunrise = 12. - hdl  !hrs, fixes time from noon till sun up by using half the day length
sundown = 12. + hdl

!find next days sunrise
hdlnext = 0.5 * dayl(2) / 3600.
sunrise_next = 12 - hdlnext

!this gets the time of peak temperature
peakt = 12. + 2. * hdl * tfpk

if (dayl(1) <= 3600.) then
  !FLAG this could be a bit of an ugly way of doing this... it could bias the weather to be mostly
  !hot since the day (which is the longest) will be set to the max temp of the day. The main reason for
  !this problem is that for gridded data we do not have the mean daily temp calced out in weathergen
  !like we do for the max and min. Will think of better ways to do this -JM Oct 29 08
      !FLAG check on this again after we fix weathergen JM Jan 10 09
  tday = tmax
  tnight = tmin

else if (dayl(1) < daysec) then !daylength is less than 24 hours

        !has a night and a day, calculate night first

          !find the length of the night
          hni = (43200. - 0.5 * dayl(1) + 43200. - 0.5 * dayl(2))  / 3600. !hrs

          b = (tfollow - t0) / sqrt(hni) !eqn 10 from Cesaraccio paper

          ti  = t0 * sundown                                            !sundown
          ti1 = t0 * (sundown + hni) + 2./3. * b * (hni**(3./2.))   !sunrise (next morn)

          tnight = (ti1 - ti) / hni

          if (dayl(1) > 0.) then  !regular night and day

                    !morning integral (ti is at sunrise, ti1 is at temperature peak time)
                    morn = sunrise - peakt

                    ti  = (tmin * sunrise) + (1. / pi) * (2. * a * morn)
                    ti1 = (tmin * peakt)   + (1. / pi) * (2. * a * morn * cos(pi/2. * (sunrise - peakt) / morn))

                    tam = (ti1 - ti) / (-morn)

                    !afternoon integral (ti is at temperature peak time, ti1 is at sundown)
                    ti  = t0 * peakt   - (1. / pi) * 8. * r * cos(pi / 8. * (-4.))
                    ti1 = t0 * sundown - (1. / pi) * 8. * r * cos(pi / 8. * (peakt - sundown - 4.))

                    tpm = (ti1 - ti) / (sundown - peakt)

                    tday = (tam + tpm) / 2.

          else

                tday = tnight      !only night, day = night

          end if

else !no night, only day

          !morning integral (ti is at sunrise, ti1 is at temperature peak time)
          morn = sunrise - peakt

          ti  = tmin * sunrise + 1. / pi * (2. * a * morn)
          ti1 = tmin * peakt   + 1. / pi * (2. * a * morn * cos(pi / 2. * (sunrise - peakt) / morn))

          tam = (ti1 - ti) / (-morn)

          !afternoon integral (t10 is at temperature peak time, ti1 is at sundown)
          ti  = t0 * peakt   - 1. / pi * 8. * r * cos(pi / 8. * (-4.))
          ti1 = t0 * sundown - 1. / pi * 8. * r * cos(pi / 8. * (peakt - sundown - 4.))


          tpm = (ti1 - ti) / (sundown - peakt)

          tday = (tam + tpm) / 2.

          tnight = tday

end if

! Calculate mean temperature from tday and tnight to stored into 'dayvars'
tmean = tday * (dayl(1) / daysec) + tnight * ((daysec - dayl(1)) / daysec)


end subroutine diurnaltemp

!-----------------------------------------------------------

subroutine humidity(grid,day)

! Calculate relative humidity from dew point temperature
! From Lawrence (2005) The relationship between relative humidity and the dewpoint temperature. A. MetSoc (https://doi.org/10.1175/BAMS-86-2-225)
! Equation (11)

use metvarsmod, only : dayvars

implicit none

integer(i4), intent(in) :: grid
integer, intent(in) :: day

real(sp), pointer :: tmean
real(sp), pointer :: tdew
real(sp), pointer :: rhum

real(sp) :: tmean_K
real(sp) :: tdew_K

real(dp), parameter :: L  = 2257000
real(dp), parameter :: Rw = 461.5

!------
! point pointers to global data
tmean => dayvars(grid,day)%tmean
tdew  => dayvars(grid,day)%tdew
rhum  => dayvars(grid,day)%rhum

tmean_K = tmean + 273.15
tdew_K  = tdew + 273.15

!------

rhum = 100. * exp(((1. - (tmean_K / tdew_K)) * (L / Rw)) / tmean_K)

if (rhum > 100.) rhum = 100.    ! When temp is lower than dew point


end subroutine humidity

!-----------------------------------------------------------

subroutine calctdew(grid,day)

use parametersmod, only : pi,d2r,r2d
use metvarsmod,    only : genvars,dayvars,gridlon,gridlat

implicit none

integer(i4), intent(in) :: grid
integer, intent(in) :: day

real(sp), parameter :: a = 1.26   ! Constant for Ep formula in Kimball et al. (1997)
real(sp), parameter :: c = 0.66   ! Constant for Ep formula in Kimball et al. (1997)

real(sp) :: yearday
real(sp) :: lon
real(sp) :: lat
real(sp) :: tmin
real(sp) :: tmax
real(sp) :: tavg
real(sp) :: prec    ! Annual rainfall in m



real(sp) :: Lv    ! Latent heat of vaporization in J K-1
real(sp) :: Pw    ! Water density in kg m-3
real(sp) :: dSVP  ! Rate of change of saturation vapour pressure in Pa K-1
real(sp) :: Ep    ! Potential evapotranspiration in kg m2 s-1
real(sp) :: EF    ! Ratio of Ep to annual precipitation (m)

real(sp) :: ampl  ! Seasonal variation in daylength in hour
real(sp) :: sunrise  ! Sunrise time
real(sp) :: sunset  ! Sunset time
real(sp) :: dayl  ! Daylength in seconds
real(sp) :: dsol  ! Solar declination angle in degrees
real(sp) :: w  ! Elevation angle in degrees
real(sp) :: Ma    ! Mean anomaly of orbit in rad
real(sp) :: va    ! True anomaly of orbit in rad
real(sp) :: Rd    ! Actual distance between sun and Earth at yearday in Gm

real(sp) :: Rn    ! Daily average insolation in W m-2
real(sp) :: Gn    ! Daily average surface conductive energy flux in W m-2

real(sp) :: Td    ! Dew point temperature in degree Celcius
real(sp) :: RH
real(sp) :: Tw

!---

tmin = dayvars(grid,day)%tmin
tmax = dayvars(grid,day)%tmax
tavg = dayvars(grid,day)%tmean
prec = sum(genvars%pre(5:16)) / 1000.   ! Annual rainfall in m

lon = gridlon(grid)
lat = gridlat(grid)
dsol = dayvars(grid,day)%dsol

yearday = day
dayl = dayvars(grid,day)%dayl * 3600.    ! daylength in seconds

!---

Lv = 1.91846 * (10**6) * (((tavg + 273.15) / (((tavg + 273.15) - 33.91))) ** 2 )    ! Henderson-Sellers (1984) in Davis et al. (2017)

Pw = 1000. * (1 - (((tavg - 3.9863) ** 2) / 508929.2) * ((tavg + 288.9414) / (tavg + 68.12963)))   ! Graf (2009)

tavg = tavg + 273.15

! dSVP = 6.1078 * exp((17.269 * tavg) / (237.3 + tavg)) * ((237.3 + tavg) * 17.269 - 17.269 * tavg) / ((237.3 + tavg) ** 2)   ! Running & Coughlan (1988)
!
! dSVP = 100. * dSVP     ! Convert mbar to Pa


tavg = (tmin + tmax) / 2.

dSVP = 0.6113 * exp(19.85356 - 5423. / tavg) * (5423. / (tavg ** 2))

dSVP = 1000 * dSVP

tavg = tavg - 273.15

!-------------------------------

! ampl = exp(7.42 + 0.045 * lat) / 3600.
!
! dayl = ampl * (sin((yearday - 79.)) * 0.01721) + 12.

!-------------------------------

Ma = (2. * pi) * (yearday - 4.) / 365.256363

va = Ma + 0.0333988 * sin(Ma) + 0.0003486 * sin(2*Ma) + 0.0000050 * sin(3*Ma)

Rd = 149.457 * (1 - 0.0167 ** 2) / (1 + 0.0167 * cos(va))

!---

dsol = 23.44 * cos(2.*pi * (yearday - 172.) / 365.)   ! Solar declination angle from Practical Meteorology Textbook

!---

dsol = dsol * d2r   ! Convert to radian

lat = lat * d2r   ! Convert to radian

lon = lon * d2r   ! Convert to radian

!---

Rn = (-1.) * tan(lat) * tan(dsol)

! Rn = min((pi / 180.), max((pi / (-180.)), Rn))

Rn = acos(Rn)

Rn = Rn * sin(lat) * sin(dsol) + cos(lat) * cos(dsol) * sin(Rn)   ! CHECKKKKKKKK
!
Rn = (1361. / pi) * ((149.457 / Rd) ** 2) * Rn

! print *, Rn, dayvars(day)%srad * 1000. / (3600. * 24)

Rn = dayvars(grid,day)%srad * 1000. / (3600. * 24)

!---

sunrise = (24. / (2 * pi)) * ((-1.) * lon + acos((sin(lat) * sin(dsol)) / (cos(lat) * cos(dsol))))

sunset = (24. / (2 * pi)) * ((-1.) * lon - acos((sin(lat) * sin(dsol)) / (cos(lat) * cos(dsol))))

if (sunrise < 0.) sunrise = sunrise + 24

if (sunset < 0.) sunset = sunset + 24

! dayl = (max(sunrise,sunset) - min(sunrise,sunset)) * 3600.

!---

Gn = 0.1 * Rn     ! Kimball et al. (1997) estimation of surface conductive energy flux

!-------------------------------

! Ep = (a * (dSVP / (dSVP + c)) * (Rn - Gn)) / Lv   ! Kimball et al. (1997)

Ep = dayvars(grid,day)%dpet / (3600. * 24.)

EF = ((Ep / Pw) * dayl) / prec     ! Kimball et al. (1997)

! print *, prec, Ep, Pw, EF

!-------------------------------
tmin = tmin + 273.15
tmax = tmax + 273.15
tavg = tavg + 273.15

Td = tmin * (-0.127 + 1.121 * (1.003 - 1.444 * EF + 12.312 * (EF**2) - 32.766 * (EF**3)) + 0.0006 * (tmax - tmin))    ! Kimball et al. (1997)

! Td = (-0.127 + 1.121 * (1.003 - 1.444 * EF + 12.312 * (EF**2) - 32.766 * (EF**3)) + 0.0006 * (tmax - tmin))    ! Kimball et al. (1997)

! Td = -0.0360 * tavg + 0.9679 * tmin + 0.0072 * (tmax - tmin) + 1.0119   ! Hubbard et al. (2003)

Td = Td - 273.15

RH = 100 * (EXP((17.625*Td)/(243.04+Td))/EXP((17.625*tavg)/(243.04+tavg)))

Tw = tavg * atan(0.151977 * sqrt(RH + 8.313659)) + atan(tavg + RH)            &
     - atan(RH - 1.676331) + 0.00391838 * (RH ** 1.5) * atan(0.023101 * RH)   &
     - 4.686035                 ! From Stull (2011) Wet-Bulb Temperature from Relative Humidity and Air Temperature


dayvars(grid,day)%tdew = Td

dayvars(grid,day)%rhum = RH

! print *, sunrise
! print *, sunset
! print *, dayl
! print *, Lv
! print *, dSVP
! print *, 'Pw:  ', Pw
! print *, 'Rn:  ', Rn
! print *, 'Ep:  ', Ep
! print *, 'EF:  ', EF * prec * 365
! print *, 'Dayl:', dayl
! print *, 'tavg:', tavg
! print *, 'tmin:', tmin - 273.15
! print *, 'DewT:', Td
! print *, 'RH:  ', RH
! print *, 'TWet:', Tw


end subroutine calctdew


end module diurnaltempmod
