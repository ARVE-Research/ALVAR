module hourlyprecmod

! This module contains subroutine to disaggreagte 24 hour total daily precipitation
! into hourly series. Disaggregation process follows basic statistical distributions
! for rough estimate of hourly values
!
! (1) Regression model used to estimate hourly rainfall statistics from daily statistics (mean, variance, skewness, dry interval)
! (2) Gamma distribution to estimate intensity of rainfall (i.e. mm h-1) at specific hour
!     ** NOTE: model currently randomly decide on wet hours instead of a statistical estimation of storm onset time
!              The simple model create similar time series pattern to more sophisticated cascade rain models at hourly scale
!              HOWEVER, this method will not be appropriate for any smaller time-step.
!
! Beuchat et al. (2011) Toward a robust method for subdaily rainflal downscaling from daily data. doi:10.1029/2010WR010342
! Menabde & Sivapalan (2000) Modeling of rainfall time series and extremes  https://doi.org/10.1029/2000WR900197

use randomdistmod, only : randomstate,ran_seed,ranu,ranur,ran_gamma

implicit none

public :: hourlyprec

type(randomstate) :: hprec_rndst

contains

!-----------------------------------------------------------------------

subroutine hourlyprec(year,grid,ndyear,day,dprec,prec,tmean,rhum,hprec)

! Subroutine written by Leo O Lai (Jul 2021)

use parametersmod, only : i4,sp,dp
use statsmod,      only : gamma_cdf_inv
use utilitiesmod,  only : getmonth

!arguments
integer(i4),               intent(in)    :: year
integer(i4),               intent(in)    :: grid
integer(i4),               intent(in)    :: ndyear
integer(i4),               intent(in)    :: day       ! Julian day (1 to 366)
real(sp),                  intent(in)    :: dprec
real(sp),    dimension(:), intent(in)    :: prec
real(sp),    dimension(:), intent(in)    :: tmean
real(sp),    dimension(:), intent(in)    :: rhum
real(dp),    dimension(:), intent(inout) :: hprec

integer(i4) :: month                 ! Current month (1 to 12)
integer(i4) :: sday                  ! Index of start day of current month
integer(i4) :: eday                  ! Index of end day of current month

! Met variables
real(sp), allocatable, dimension(:) :: mprec     ! Monthly array of daily 24h total precipitation (mm)
real(sp), allocatable, dimension(:) :: mtmean    ! Monthly array of daily mean temperature (C)
real(sp), allocatable, dimension(:) :: mrhum     ! Monthly array of daily mean humidity (%)

! Hourly precipitation regression variables for gamma random number generator
real(sp) :: mean_1             ! Mean of hourly precipitation
real(sp) :: var_1              ! Variance of hourly precipitation
real(sp) :: skew_1             ! Skewness of hourly precipitation
real(sp) :: dryi_1             ! Dry interval of hourly precipitation

integer(i4) :: dry_hour        ! Number of dry hours in current day
integer(i4) :: wet_hour        ! Number of wet hours in current day

real(dp) :: shape              ! Shape parameter of gamma distribution
real(dp) :: scale              ! Scale parameter of gamma distribution

! Other variables
real(sp),    dimension(24) :: ran
real(sp),    dimension(24) :: sort_val
integer(i4), dimension(24) :: sort_loc
real(sp)                   :: val
integer(i4)                :: loc

real(sp) :: diff
integer :: i

!------
! Get daily met variables of the current month, calculate start and end array index of the month

call getmonth(day,ndyear,month,sday,eday)

allocate(mprec(eday-sday+1))
allocate(mtmean(eday-sday+1))
allocate(mrhum(eday-sday+1))

mprec  = prec(sday:eday)
mtmean = tmean(sday:eday)
mrhum  = rhum(sday:eday)

!------

! Initialize the random state for hourly prec simulation
if (year == 1 .and. grid == 1 .and. day == 1) call ran_seed(0, hprec_rndst)

!------

if (dprec == 0.) then      ! Dry day

  hprec = 0.

else                       ! Wet day

  ! Initialize hprec = 0 due to dry interval / dry hours
  hprec = 0.

  ! Get hourly rainfall statistics (variance, skewness and wet interval) using disaggregation
  ! regression model from Beuchat et al. (2011)
  if (dprec == sum(mprec)) then    ! only one wet day in month

    var_1  = 0.1
    skew_1 = 0.01
    dryi_1 = 0.5

  else

    call hprec_regression(mprec,mtmean,mrhum,var_1,skew_1,dryi_1)

  end if

  !------

  ! Number of wet hours in a day
  dry_hour = floor(dryi_1 * 24.)

  wet_hour = 24 - dry_hour

  ! Mean hourly rainfall across all wet hours
  mean_1 = dprec / wet_hour

  ! Shape and scale parameters of gamma function for estimating hourly rainfall depth
  shape = dble((mean_1 ** 2) / var_1)
  scale = dble(mean_1 / var_1)

  !------

  ! Generate 24 random numbers between 0 and 1 to decide on wet hours
  do i = 1, 24
    ran(i) = ranur(hprec_rndst)
  end do

  !------
  ! Sort the random [0,1] values in decending order to randomly distribute wet hours across the day
  ! E.g. if wet_hour = 10, the highest 10 values in ran(i) will be taken as the 10 hours where
  ! rainfall is observed, while the rest of the 14 hours are assumed to be dry
  !
  ! The values are then fitted into a gamma function to generate a random gamma number as a measure of
  ! hourly rain depth
  call sort_descend(ran,sort_val,sort_loc)

  do i = 1, wet_hour

    val = sort_val(i)
    loc = sort_loc(i)

    !---

    hprec(loc) = ran_gamma_num(shape,scale,hprec_rndst)

    ! Round to zero if hourly prec is less than 1.e-2
    if (hprec(loc) < 0.01) hprec(loc) = 0.

  end do

  !------

  ! if (dprec > 20.) print *, sum(hprec), dprec, dprec - sum(hprec), wet_hour

  ! Compare and check hourly precipitation with daily 24h total value
  diff = sum(hprec) / dprec

  ! Compensate for difference created by random gamma number generator
  where (hprec > 0.)
    hprec = hprec / diff
  end where

end if

deallocate(mprec)
deallocate(mtmean)
deallocate(mrhum)

end subroutine hourlyprec


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


subroutine hprec_regression(prec,tmean,rhum,var_1,skew_1,dryi_1)

! Get hourly rainfall statistics (variance, skewness and wet interval) using disaggregation
! regression model from Beuchat et al. (2011)

use parametersmod, only : i4,sp

implicit none

real(sp), parameter                :: intcpt_var  = 0.0384
real(sp), parameter, dimension(16) :: coeff_var   = [0.805, -0.878, 0.0561, -0.0219, &
                                                    -1.03, -0.187, 0.308, -0.000886, &
                                                    -0.0137, -0.000101, -0.000389, 0.016, &
                                                    0.035, 0.00207, -0.00551, -0.00548]

real(sp), parameter                :: intcpt_weti = 1.62
real(sp), parameter, dimension(17) :: coeff_weti  = [0.871, -0.749, 0.0212, -0.0165, &
                                                    -0.538, 0.478, -0.394, 0.325, &
                                                    0.0169, -0.252, 0.290, 0.0000932, &
                                                    0.00272, -0.00141, -0.000000471, -0.00000204, &
                                                    0.192]

real(sp), parameter                :: intcpt_skew = 2.16
real(sp), parameter, dimension(16) :: coeff_skew  = [0.615, -0.289, 0.0263, -0.435, &
                                                    0.425, -0.0199, -0.917, 0.100, &
                                                    -0.188, 0.00645, -0.0848, -0.0629, &
                                                    -0.00000569, -0.0000225, 0.00934, 0.000656]

real(sp), dimension(:), intent(in)  :: prec
real(sp), dimension(:), intent(in)  :: tmean
real(sp), dimension(:), intent(in)  :: rhum
real(sp),               intent(out) :: var_1
real(sp),               intent(out) :: skew_1
real(sp),               intent(out) :: dryi_1

real(sp) :: mean_24
real(sp) :: var_24
real(sp) :: sd_24
real(sp) :: skew_24
real(sp) :: auto_24
real(sp) :: dryi_24
real(sp) :: dryi_48
real(sp) :: elev
real(sp) :: tas
real(sp) :: hur
real(sp) :: vartas_24

real(sp), dimension(:), allocatable :: terms

integer(i4) :: nd
integer(i4) :: dry_1day
integer(i4) :: dry_2day
integer(i4) :: i

!------

nd = size(prec)

!------

dry_1day = 0
dry_2day = 0

do i = 1, nd
  if (prec(i) == 0.) dry_1day = dry_1day + 1
end do

do i = 1, (nd-1)
  if (prec(i) == 0. .and. prec(i+1) == 0.) dry_2day = dry_2day + 1
end do

!------

mean_24   = sum(prec) / nd
var_24    = variance(prec)
sd_24     = sqrt(var_24)
skew_24   = skewness(prec,sd_24)
auto_24   = autocorrel(prec,1)
dryi_24   = dry_1day / nd
dryi_48   = dry_2day / (nd-1)
elev      = 100.
tas       = sum(tmean) / nd
hur       = sum(rhum) / nd
vartas_24 = variance(tmean)

! Log-transform variance and skewness for model (Beuchat et al., 2011)
mean_24 = log(mean_24)
var_24  = log(var_24)
skew_24 = log(abs(skew_24))       ! Temporary solution to filter occurance of (very small) negative values of skewness

!-----------------------------------------------------------------------
! Hourly variance from 24 hours statistics regression model
allocate(terms(16))

terms(1)  = var_24 - 5.09
terms(2)  = 5.09 - var_24
terms(3)  = tas - 3.21
terms(4)  = 3.21 - tas
terms(5)  = auto_24 - 0.103
terms(6)  = skew_24 - 1.67
terms(7)  = 1.67 - skew_24
terms(8)  = (tas - 3.21) * (hur - 66.3)
terms(9)  = (dryi_48 - 0.141) * (tas - 3.21)
terms(10) = elev - 191.1
terms(11) = 191.1 - elev
terms(12) = (5.09 - var_24) * (tas - 17.3)
terms(13) = (5.09 - var_24) * (skew_24 - 2.60) * (17.4 - tas)
terms(14) = (5.09 - var_24) * (2.61 - skew_24) * (17.4 - tas)
terms(15) = hur - 48.4
terms(16) = 48.4 - hur

!---

if (tas - 3.21 < 0) then

  terms(8) = 0
  terms(9) = 0

else if (hur - 66.3 < 0) then

  terms(8) = 0

else if (dryi_48 - 0.141 < 0) then

  terms(9) = 0

else if (5.09 - var_24 < 0) then

  terms(12) = 0
  terms(13) = 0
  terms(14) = 0

else if (17.4 - tas < 0) then

  terms(13) = 0
  terms(14) = 0

else if (skew_24 - 2.60 < 0) then

  terms(13) = 0

else if (2.61 - skew_24 < 0) then

  terms(14) = 0

end if

!---

do i = 1, 16
  if (terms(i) < 0.) terms(i) = 0.
end do

!---
! Linear regression sum of intercept and all terms with log-transformation
var_1 = intcpt_var + sum(coeff_var * terms)

var_1 = exp(var_1)

deallocate(terms)

!-----------------------------------------------------------------------
! Hourly wet interval from 24 hours statistics regression model
allocate(terms(17))

terms(1)  = dryi_24 - 0.0297
terms(2)  = 0.0297 - dryi_24
terms(3)  = tas - 1.85
terms(4)  = 1.85 - tas
terms(5)  = mean_24 - 1.37
terms(6)  = 1.37 - mean_24
terms(7)  = dryi_48 - 0.187
terms(8)  = -0.187 - dryi_48
terms(9)  = tas - 16.7
terms(10) = auto_24 - 0.145
terms(11) = 0.145 - auto_24
terms(12) = (tas - 1.85) * (65.5 - hur)
terms(13) = (var_24 - 1.26) * (tas - 1.85)
terms(14) = (1.26 - var_24) * (tas - 1.85)
terms(15) = (elev - 210.9) * (tas - 1.85) * (hur - 65.5)
terms(16) = (210.9 - elev) * (tas - 1.85) * (hur - 65.5)
terms(17) = (mean_24 - 1.087) * (dryi_24 - 0.0296)

!---

if (tas - 1.85 < 0) then

  terms(12) = 0
  terms(13) = 0
  terms(14) = 0
  terms(15) = 0
  terms(16) = 0

else if (hur - 65.5 < 0) then

  terms(12) = 0
  terms(15) = 0
  terms(16) = 0

else if (var_24 - 1.26 < 0) then

  terms(13) = 0

else if (1.26 - var_24 < 0) then

  terms(14) = 0

else if (elev - 210.9 < 0) then

  terms(15) = 0

else if (210.9 - elev < 0) then

  terms(16) = 0

else if (mean_24 - 1.087 < 0) then

  terms(17) = 0

else if (dryi_24 - 0.0296 < 0) then

  terms(17) = 0

end if

!---

do i = 1, 17
  if (terms(i) < 0.) terms(i) = 0.
end do

!---
! Linear regression sum of intercept and all terms
! logit-transform output weti (Beuchat et al., 2011)
dryi_1 = intcpt_weti + sum(coeff_weti * terms)
dryi_1 = exp(dryi_1) / (1 + exp(dryi_1))

deallocate(terms)

!-----------------------------------------------------------------------
! Hourly skewness from 24 hours statistics regression model
allocate(terms(16))

terms(1)  = skew_24 - 1.52
terms(2)  = 1.52 - skew_24
terms(3)  = tas - 0.0345
terms(4)  = mean_24 - 1.17
terms(5)  = 1.17 - mean_24
terms(6)  = hur - 70.2
terms(7)  = auto_24 - 0.114
terms(8)  = var_24 - 3.05
terms(9)  = 3.05 - var_24
terms(10) = (tas - 11.0) * (hur - 70.2)
terms(11) = (skew_24 - 1.52) * (vartas_24 - 6.44)
terms(12) = (skew_24 - 1.52) * (6.44 - vartas_24)
terms(13) = (elev - 203.9) * (tas - 11) * (hur - 70.2)
terms(14) = (203.9 - elev) * (tas - 11) * (hur - 70.2)
terms(15) = (skew_24 - 2.35) * (tas - 11.0) * (hur - 70.2)
terms(16) = (2.35 - skew_24) * (tas - 11.0) * (hur - 70.2)

!---

if (tas - 11.0 < 0) then

  terms(10) = 0
  terms(13) = 0
  terms(14) = 0
  terms(15) = 0
  terms(16) = 0

else if (hur - 70.2 < 0) then

  terms(10) = 0
  terms(13) = 0
  terms(14) = 0
  terms(15) = 0
  terms(16) = 0

else if (skew_24 - 1.52 < 0) then

  terms(11) = 0
  terms(12) = 0

else if (vartas_24 - 6.44 < 0) then

  terms(11) = 0

else if (6.44 - vartas_24 < 0) then

  terms(12) = 0

else if (elev - 203.9 < 0) then

  terms(13) = 0

else if (203.9 - elev < 0) then

  terms(14) = 0

else if (skew_24 - 2.35 < 0) then

  terms(15) = 0

else if (2.35 - skew_24 < 0) then

  terms(16) = 0

end if

!---

do i = 1, 16
  if (terms(i) < 0.) terms(i) = 0.
end do

!---
! Linear regression sum of intercept and all terms
! logit-transform output weti (Beuchat et al., 2011)
skew_1 = intcpt_skew + sum(coeff_skew * terms)
skew_1 = exp(skew_1)

deallocate(terms)

end subroutine hprec_regression

!-----------------------------------------------------------------------

real(sp) function variance(a)

! Find variance of a real array

use parametersmod, only : i4,sp

implicit none

real(sp), dimension(:), intent(in) :: a

real(sp) :: summ
real(sp) :: mean
integer(i4) :: len
integer(i4) :: i

!---

len = size(a)
mean = sum(a) / len

summ = 0.

do i = 1, len

  summ = summ + (a(i) - mean) ** 2.

end do

variance = summ / len

end function variance

!-----------------------------------------------------------------------

real(sp) function skewness(a,sd)

! Find the skewness of a real array with given standard deviation (sd)

use parametersmod, only : i4,sp

implicit none

real(sp), dimension(:), intent(in) :: a
real(sp)              , intent(in) :: sd

real(sp) :: summ
real(sp) :: mean
integer(i4) :: len
integer(i4) :: i

!---

len = size(a)
mean = sum(a) / len

summ = 0.

do i = 1, len

  summ = summ + ((a(i) - mean) ** 3.)

end do

skewness = summ / ((real(len) - 1.) * (sd ** 3.))

end function skewness

!-----------------------------------------------------------------------

real(sp) function autocorrel(a,k)

! Find the auto-correlation of real array with given lag-k

use parametersmod, only : i4,sp

implicit none

real(sp), dimension(:), intent(in) :: a
integer               , intent(in) :: k

real(sp) :: summ
real(sp) :: summ2
real(sp) :: mean
integer(i4) :: len
integer(i4) :: i

!---

len = size(a)
mean = sum(a) / len

summ = 0.

do i = 1, (len-k)

  summ = summ + ((a(i) - mean) * (a(i+k) - mean))

end do

!---

summ2 = 0.

do i = 1, len

  summ2 = summ2 + (a(i) - mean) ** 2.

end do

!---

autocorrel = summ / summ2

end function autocorrel

!-----------------------------------------------------------------------

subroutine sort_ascend(a,sort_val,sort_loc)

! A subroutine to sort a real array in decending order
! Return two arrays for decending values and the index locations

use parametersmod, only : i4,sp

implicit none

real(sp),    dimension(:), intent(in)  :: a
real(sp),    dimension(:), intent(out) :: sort_val
integer(i4), dimension(:), intent(out) :: sort_loc

real(sp), dimension(:), allocatable :: b

real(sp)    :: small
real(sp)    :: val
integer(i4) :: loc
integer(i4) :: i

!------

! Create a copy of the input array for manipulation
allocate(b(size(a)))

b = a

small = minval(b)

!------

do i = 1, size(b)

  val = maxval(b)
  loc = maxloc(b, dim=1)

  sort_val(i) = val
  sort_loc(i) = loc

  b(loc) = small - 1.

end do

end subroutine sort_ascend

!-----------------------------------------------------------------------

subroutine sort_descend(a,sort_val,sort_loc)

! A subroutine to sort a real array in decending order
! Return two arrays for decending values and the index locations

use parametersmod, only : i4,sp

implicit none

real(sp),    dimension(:), intent(in)  :: a
real(sp),    dimension(:), intent(out) :: sort_val
integer(i4), dimension(:), intent(out) :: sort_loc

real(sp), dimension(:), allocatable :: b

real(sp)    :: large
real(sp)    :: val
integer(i4) :: loc
integer(i4) :: i

!------

! Create a copy of the input array for manipulation
allocate(b(size(a)))

b = a

large = maxval(b)

!------

do i = 1, size(b)

  val = minval(b)
  loc = minloc(b, dim=1)

  sort_val(i) = val
  sort_loc(i) = loc

  b(loc) = large + 1.

end do

end subroutine sort_descend

!-----------------------------------------------------------------------

FUNCTION ran_gamma_num(G,H, rndst)
! *     Function copied from https://www.ucl.ac.uk/~ucakarc/work/software/randgen.f (NEED REFERENCING!!!!!!)
! *
! *       Returns a random number with a gamma distribution with mean
! *       G/H and variance G/(H^2). (ie. shape parameter G & scale
! *       parameter H)
! *

use parametersmod, only : dp

real(dp) :: C,D,R,ran_gamma_num,ZBQLU01,G,H,A,z1,z2,B1,B2,M
real(dp) :: U1,U2,U,V,TEST,X
real(dp) :: c1,c2,c3,c4,c5,w

type(randomstate), intent(inout) :: rndst

ran_gamma_num = 0.0D0

IF ( (G.LE.0.0D0).OR.(H.LT.0.0D0) ) THEN
 WRITE(*,*) "ERMMMMM"
 RETURN
ENDIF

IF (G.LT.1.0D0) THEN
889    u=ranur(rndst)
 v=ranur(rndst)
 if (u.gt.exp(1.0d0)/(g+exp(1.0d0))) goto 891
 g = max(g,1e-2_dp)
 u = max(u,1e-2_dp)
 ran_gamma_num=((g+exp(1.0d0))*u/exp(1.0d0))**(1.0d0/g)
 if (v.gt.exp(-ran_gamma_num)) then
goto 889
 else
goto 892
 endif
891    ran_gamma_num=-log((g+exp(1.0d0))*(1.0d0-u)/(g*exp(1.0d0)))
 if (v.gt.ran_gamma_num**(g-1.0)) goto 889
892    ran_gamma_num=ran_gamma_num/h
 RETURN
ELSEIF (G.LT.2.0D0) THEN
 M = 0.0D0
elseif (g.gt.10.0d0) then
 c1=g-1.0d0
 c2=(g-1.0d0/(6.0d0*g))/c1
 c3=2.0d0/c1
 c4=c3+2.0d0
 c5=1.0d0/sqrt(g)
777    u=ranur(rndst)
 v=ranur(rndst)
 if (g.gt.2.50d0) then
u=v+c5*(1.0d0-1.860d0*u)
 endif
 if (u.le.0.0d0.or.u.ge.1.0d0) goto 777
 w=c2*v/u
 if (c3*u+w+1.0d0/w.le.c4) goto 778
 if (c3*log(u)-log(w)+w.ge.1.0d0) goto 777
778    ran_gamma_num=c1*w/h
 return
ELSE
 M = -(G-2.0D0)
ENDIF
R = 0.50D0
a = ((g-1.0d0)/exp(1.0d0))**((g-1.0d0)/(r+1.0d0))
C = (R*(M+G)+1.0D0)/(2.0D0*R)
D = M*(R+1.0D0)/R
z1 = C-DSQRT(C*C-D)
! *
! *     On some systems (e.g. g77 0.5.24 on Linux 2.4.24), C-DSQRT(C*C)
! *     is not exactly zero - this needs trapping if negative.
! *
IF ((Z1-M.LT.0.0D0).AND.(Z1-M.GT.-1.0D-12)) Z1 = M
z2 = C+DSQRT(C*C-D)
B1=(z1*(z1-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z1-M)/(R+1.0D0))
B2=(z2*(z2-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z2-M)/(R+1.0D0))
50    U1=ranur(rndst)
U2=ranur(rndst)
U=A*U1
V=B1+(B2-B1)*U2
X=V/(U**R)
IF (X.LE.M) GOTO 50
! x = max(x,1e-2_dp)
! m = max(m,1e-2_dp)
! g = max(g,1e-2_dp)
! r = max(r,1e-2_dp)
TEST = ((X-M)**((G-1)/(R+1)))*EXP(-(X-M)/(R+1.0D0))
IF (U.LE.TEST) THEN
 ran_gamma_num = (X-M)/H
ELSE
 GOTO 50
ENDIF
!   FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
! +' ran_gamma_num',/5X, '(both parameters must be positive)',/)

END
! ***************************************************************


end module hourlyprecmod
