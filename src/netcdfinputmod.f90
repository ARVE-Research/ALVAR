module netcdfinputmod

use parametersmod, only : i2,i4,sp,dp
use metvarsmod,    only : xlen,ylen,ilen,tlen,lon,lat,indx,time,srt,cnt,cntt,p0,p1,monvars
use simplesoilmod, only : nl,zpos,soilvars
use outputmod,     only : infompi
use errormod,      only : ncstat,netcdf_err
use getdatamod,    only : readdata
use netcdf
use mpi

implicit none

character(100), parameter :: soilfile = './soil30m_list-formatted.nc'

contains

!---------------------------------------------------------------------

subroutine metdatainput(info)

! Read in the full array (all months for all gridcells) of monthly variables from netcdf file

type(infompi), target, intent(in) :: info

character(100), pointer :: infile
integer :: ifid
integer :: dimid
integer :: varid

integer(i4) :: gridcount
integer(i4), dimension(2) :: start
integer(i4), dimension(2) :: count

! monthly input driver variables
real(sp), allocatable, dimension(:,:) :: tmp        ! mean monthly temperature (degC)
real(sp), allocatable, dimension(:,:) :: dtr        ! mean monthly diurnal temperature range (degC)
real(sp), allocatable, dimension(:,:) :: pre        ! total monthly precipitation (mm)
real(sp), allocatable, dimension(:,:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
real(sp), allocatable, dimension(:,:) :: cld        ! mean monthly cloud cover (fraction)
real(sp), allocatable, dimension(:,:) :: wnd        ! mean monthly 10m windspeed (m s-1)


integer :: i


!----------------------------------------------------
! Pointers to mpi derived type variables
infile => info%infile

gridcount = cnt(1)

! -----------------------------------------------------
! INPUT: Read dimension IDs and lengths of dimensions

ncstat = nf90_open(infile,nf90_nowrite,ifid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           ! Open netCDF-file (inpput file name, no writing rights, assigned file number)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! Check for errors (after every step)

ncstat = nf90_inq_dimid(ifid,'lon',dimid)              ! get dimension ID from dimension 'lon' in the input file
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lat',dimid)              ! get dimension ID from dimension 'lon' in the input file
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'index',dimid)             ! get dimension ID from dimension 'index' in the input file
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ilen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'time',dimid)             ! Get dimension ID for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)   ! Get length of dimension 'time' and assign it to variable tlen
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----------------------------------------------------
! Read variable IDs and values

allocate(lon(xlen))       ! Allocate length to longitude array
allocate(lat(ylen))       ! Allocate length to latitude array
allocate(time(tlen))      ! Allocate length to time array
allocate(indx(xlen,ylen)) ! Allocate length to index array

ncstat = nf90_inq_varid(ifid,"lon",varid)                ! Get variable ID for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lon)                    ! Get variable values for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)                ! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lat)                    ! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"index",varid)               ! Get variable ID for index
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,indx)                    ! Get variable values for index
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"time",varid)               ! Get variable ID for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,time)                   ! Get variable values for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----------------------------------------------------
! Read in variables from netcdf file by each gridcell
srt = srt
cnt = cnt

allocate(tmp(cnt(1),cntt))
allocate(dtr(cnt(1),cntt))
allocate(pre(cnt(1),cntt))
allocate(wet(cnt(1),cntt))
allocate(cld(cnt(1),cntt))
allocate(wnd(cnt(1),cntt))

call readdata(ifid,'tmp',srt,cnt,tmp(:,p0:p1))
call readdata(ifid,'dtr',srt,cnt,dtr(:,p0:p1))
call readdata(ifid,'pre',srt,cnt,pre(:,p0:p1))
call readdata(ifid,'wet',srt,cnt,wet(:,p0:p1))
call readdata(ifid,'cld',srt,cnt,cld(:,p0:p1))
call readdata(ifid,'wnd',srt,cnt,wnd(:,p0:p1))

!---------------------------------------------------------------------
! copy first and last year into buffers at each end
if (p0 == 13) then
  tmp(:,1:12) = tmp(:,13:24)
  dtr(:,1:12) = dtr(:,13:24)
  pre(:,1:12) = pre(:,13:24)
  wet(:,1:12) = wet(:,13:24)
  cld(:,1:12) = cld(:,13:24)
  wnd(:,1:12) = wnd(:,13:24)
end if

if (cntt > p1) then
  tmp(:,p1+1:cntt) = tmp(:,p1-11:p1)
  dtr(:,p1+1:cntt) = dtr(:,p1-11:p1)
  pre(:,p1+1:cntt) = pre(:,p1-11:p1)
  wet(:,p1+1:cntt) = wet(:,p1-11:p1)
  cld(:,p1+1:cntt) = cld(:,p1-11:p1)
  wnd(:,p1+1:cntt) = wnd(:,p1-11:p1)
end if

!---------------------------------------------------------------------
! enforce reasonable values of prec and wetdays

where (pre > 0.1)
  wet = max(wet,1.)
  wet = min(wet,10. * pre)
elsewhere
  pre = 0.
  wet = 0.
end where

!---------------------------------------------------------------------
! Copy monthly time series into the mon_metvar derived type structure
do i = 1, cnt(1)

  monvars(i)%tmp(:) = tmp(i,:)
  monvars(i)%dtr(:) = dtr(i,:)
  monvars(i)%pre(:) = pre(i,:)
  monvars(i)%wet(:) = wet(i,:)
  monvars(i)%cld(:) = cld(i,:)
  monvars(i)%wnd(:) = wnd(i,:)

end do

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


end subroutine metdatainput

!---------------------------------------------------------------------

subroutine soildatainput(info)

! Read in the full array all gridcells and all layers, of soil variales (sand, clay and cfvo composition)

type(infompi), target, intent(in) :: info

integer :: sfid
integer :: dimid
integer :: varid

integer(i4) :: gridstart
integer(i4) :: gridcount

integer(i2), allocatable, dimension(:,:) :: var_in
real(sp),    allocatable, dimension(:,:) :: values

real(sp)    :: scale_factor
real(sp)    :: add_offset
integer(i2) :: missing_value

integer :: i

!------

gridstart = srt(1)
gridcount = cnt(1)

! -----------------------------------------------------
! INPUT: Read in layer depth variable (zpos)

ncstat = nf90_open(soilfile,nf90_nowrite,sfid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           !
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(sfid,"zpos",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(sfid,varid,zpos)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! -----------------------------------------------------
! INPUT: Read in sand composition

allocate(var_in(gridcount,nl))
allocate(values(gridcount,nl))

ncstat = nf90_inq_varid(sfid,"sand",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(sfid,varid,var_in,start=[gridstart,1],count=[gridcount,nl])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(sfid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(sfid,varid,"scale_factor",scale_factor)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value)
  values = real(var_in) * scale_factor
elsewhere
  values = -9999.
end where

do i = 1, gridcount
  soilvars(i)%sand(1:nl) = values(i,1:nl)
end do

deallocate(var_in)
deallocate(values)

! -----------------------------------------------------
! INPUT: Read in sand composition

allocate(var_in(gridcount,nl))
allocate(values(gridcount,nl))

ncstat = nf90_inq_varid(sfid,"clay",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(sfid,varid,var_in,start=[gridstart,1],count=[gridcount,nl])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(sfid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(sfid,varid,"scale_factor",scale_factor)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value)
  values = real(var_in) * scale_factor
elsewhere
  values = -9999.
end where

do i = 1, gridcount
  soilvars(i)%clay(1:nl) = values(i,1:nl)
end do

deallocate(var_in)
deallocate(values)

! -----------------------------------------------------
! INPUT: Read in sand composition

allocate(var_in(gridcount,nl))
allocate(values(gridcount,nl))

ncstat = nf90_inq_varid(sfid,"cfvo",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(sfid,varid,var_in,start=[gridstart,1],count=[gridcount,nl])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(sfid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(sfid,varid,"scale_factor",scale_factor)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value)
  values = real(var_in) * scale_factor
elsewhere
  values = -9999.
end where

do i = 1, gridcount
  soilvars(i)%cfvo(1:nl) = values(i,1:nl)
end do

deallocate(var_in)
deallocate(values)

!------

ncstat = nf90_close(sfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


end subroutine soildatainput


end module netcdfinputmod
