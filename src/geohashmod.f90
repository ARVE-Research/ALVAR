module geohashmod

! Module for setting the seed based upon the location on the earth

use parametersmod, only : sp,dp,i4

implicit none

integer(i4), parameter :: mv = huge(i4)  !add offset to fit into values into signed number

contains

!------

integer(i4) function geohash(lon,lat)

! encode a lon-lat pair into a signed 4-byte integer

! original use: to seed a random number generator in a predictable but geographically non-uniform way
! this should be precise and provide a unique number down to 30 arc-second resolution (1km)
! Jed O. Kaplan, 2011

implicit none

! arguments

real(dp), intent(in) :: lon ! the longitude
real(dp), intent(in) :: lat ! the latitude

! parameters

real(sp),    parameter :: scale  = 120.d0      !scale factor for geohash, the larger the number the more unique values
! real(dp),    parameter :: scale  = 7200.d0      !scale factor for geohash, larger number for very high resolution dataset
real(sp),    parameter :: offset =   0.5d0     !offset to calculate pixel number assuming gridcell center coordinates are given
integer(i4), parameter :: rowlen = nint(scale * 360.d0)

! local variables

integer(i4) :: i
integer(i4) :: j

!---

i = nint(offset + scale * (lon + 180.d0))
j = nint(offset + scale * (lat + 90.d0))

geohash = i + rowlen * (j-1) - mv

end function geohash

!------

end module geohashmod
