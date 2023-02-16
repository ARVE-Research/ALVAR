module printgridmod

implicit none

contains

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

subroutine readgrdlist(gridlistfile,plon,plat)

use parametersmod, only : i4,sp
use coordsmod,     only : index,parsecoords

implicit none

character(100),                            intent(in)  :: gridlistfile
real(sp),       allocatable, dimension(:), intent(out) :: plon
real(sp),       allocatable, dimension(:), intent(out) :: plat

! Local variables
character(50), allocatable, dimension(:) :: gridlonlat
character(50) :: tmp
type(index)   :: timevals
integer(i4)   :: i
integer(i4)   :: n
integer(i4)   :: stat

!---------------------

open(10,file=gridlistfile)

!---------------------
! Count the number of lines in the gridlist file

n = 0
stat = 0

do
  read(10,*,iostat=stat) tmp
  if (stat /= 0) exit
  n = n + 1
end do

!---------------------
! Allocate variables and rewind the opened file

allocate(gridlonlat(n))
allocate(plon(n))
allocate(plat(n))

rewind(10)

!---------------------

do i = 1, n

  read(10,'(A)') gridlonlat(i)

  call parsecoords(gridlonlat(i),timevals)

  plon(i) = real(timevals%minlon)
  plat(i) = real(timevals%minlat)

  if (plon(i) > 180. .OR. plon(i) < -180) then
    write(0,*) ' '
    write(0,*) 'User-specified longitude exceeded dimension in line ', i, ' of ', gridlistfile
    write(0,*) ' '
    stop
  end if

  if (plat(i) > 90. .OR. plat(i) < -90) then
    write(0,*) ' '
    write(0,*) 'User-specified latitude exceeded dimension in line', i, ' of ', gridlistfile
    write(0,*) ' '
    stop
  end if

end do

end subroutine readgrdlist

!---------------------------------------------------------------------

subroutine setgrid(lon,lat,plon,plat,ngrid,gruns)

use parametersmod, only : i4,sp,dp
use utilitiesmod,  only : which

implicit none

real(dp),                 dimension(:), intent(in)  :: lon
real(dp),                 dimension(:), intent(in)  :: lat
real(sp),                 dimension(:), intent(in)  :: plon
real(sp),                 dimension(:), intent(in)  :: plat
integer(i4),                            intent(out) :: ngrid
integer(i4), allocatable, dimension(:), intent(out) :: gruns

! Local variables for finding grid number via location of the lon/lat arrays
real(sp),    allocatable, dimension(:) :: lon_sp
integer(i4), allocatable, dimension(:) :: lon_match
integer(i4) :: ll
integer(i4) :: i
integer(i4) :: j

!---------------------

allocate(lon_sp(size(lon)))

lon_sp = sngl(lon)

!---------------------

ngrid = 0

do i = 1, size(plon)

  call which(lon_sp,plon(i),lon_match)

  !------------------

  do j = 1, size(lon_match)

    ll = lon_match(j)

    if (lat(ll) == plat(i)) then

      ngrid = ngrid + 1

    end if

  end do

  !------------------

  deallocate(lon_match)

end do

!---------------------
! Second iteration to save the idx value if match into gruns
allocate(gruns(ngrid))

ngrid = 0

do i = 1, size(plon)

  lon_sp = sngl(lon)

  call which(lon_sp,plon(i),lon_match)

  !------------------

  do j = 1, size(lon_match)

    ll = lon_match(j)

    if (lat(ll) == plat(i)) then

      ngrid = ngrid + 1

      gruns(ngrid) = ll

    end if

  end do

  !------------------

  deallocate(lon_match)

end do

end subroutine setgrid

!---------------------------------------------------------------------

subroutine gentextfile(lon,lat,gruns,textfilenames)

use parametersmod, only : i4,sp,dp

implicit none

real(dp),    dimension(:), intent(in)    :: lon
real(dp),    dimension(:), intent(in)    :: lat
integer(i4), dimension(:), intent(in)    :: gruns

character(50), allocatable, dimension(:),intent(inout) :: textfilenames

integer(i4) :: g
integer(i4) :: n
integer(i4) :: i
integer(i4) :: j

!---------------------

allocate(textfilenames(size(gruns)))

do g = 1, size(gruns)

  n = gruns(g)

  i = int(lon(n) * 100.)
  j = int(lat(n) * 100.)

  if (i < 0) i = abs(i - 900000)
  if (j < 0) j = abs(j - 900000)

  write(textfilenames(g), "('output_',I6.6,'-',I6.6,'.txt')") i,j

  open (110, file=textfilenames(g), status='new')

  close(110)

end do

end subroutine gentextfile

!---------------------------------------------------------------------

subroutine printgrid(calcyrs,year,day,i,grid,textfilename,sv)

use parametersmod, only : i4,sp,Tfreeze
use statevarsmod,  only : statevars
use weathergenmod, only : roundto

implicit none

integer(i4),     intent(in) :: calcyrs
integer(i4),     intent(in) :: year
integer(i4),     intent(in) :: day
integer(i4),     intent(in) :: i
integer(i4),     intent(in) :: grid
character(50),   intent(in) :: textfilename
type(statevars), intent(in) :: sv

real(sp) :: soilm
real(sp) :: z
real(sp) :: dn
real(sp) :: prec
integer(i4) :: nl

!---------------------

open(110, file=textfilename, status='old', position='append')

! write(110,*) year, day, i, &
!              sv%dayvars%prec(day), &
!              sv%soilvars%Tliq / sv%soilvars%Tsat, &
!              sv%soilvars%Tsoil - 273.15

do nl = 0, 7

  if (nl == 0) then

    soilm = sv%soilvars%Tliq(1) / sv%soilvars%Tsat(1)
    ! soilm = sv%soilvars%swf(1)
    soilm = sv%soilvars%Tsoil(1) - Tfreeze

    soilm = roundto(soilm, 3)
    z = roundto(sv%soilvars%zipos(1), 3)

  else if (nl == 7) then

    soilm = sv%soilvars%Tliq(6) / sv%soilvars%Tsat(6)
    ! soilm = sv%soilvars%swf(6)
    soilm = sv%soilvars%Tsoil(6) - Tfreeze

    soilm = roundto(soilm, 3)
    z = roundto(sv%soilvars%zipos(7), 3)

  else

    soilm = sv%soilvars%Tliq(nl) / sv%soilvars%Tsat(nl)
    ! soilm = sv%soilvars%swf(nl)
    soilm = sv%soilvars%Tsoil(nl) - Tfreeze

    soilm = roundto(soilm, 3)
    z = roundto(sv%soilvars%zpos(nl), 3)

  end if

  if (i == 1) then
    dn = real(day)
    prec = sv%dayvars%dayprec(day)
  else
    dn = real(day) + 0.5
    prec = sv%dayvars%nightprec(day)
  end if

  if (year == calcyrs .and. i == 2) write(110,*) year, dn, &
                                          -z, &
                                          soilm, &
                                          prec, &
                                          sv%dayvars%tmean(day)
                                          ! sv%soilvars%zsno

  110 format (I1,1x,I3,1x,F6.3,1x,F5.3)

end do

close(110)

end subroutine printgrid

end module printgridmod
