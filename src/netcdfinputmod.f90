module netcdfinputmod

use parametersmod, only : i2,i4,sp,dp
use errormod,      only : ncstat,netcdf_err
use getdatamod,    only : readdata
use netcdf
use mpi

implicit none


contains

!---------------------------------------------------------------------

subroutine metdatainput(infile,gridcount,srt,cnt,cntt,p0,p1,xlen,ylen,tlen, &
                        inlon,inlat,time,indx,elev,invars)

! Read in the full array (all months for all gridcells) of monthly variables from netcdf file
use parametersmod, only : i2,i4,sp,dp
use statevarsmod,  only : in_metvars
use getdatamod,    only : readdata

implicit none

character(200),                                intent(in)    :: infile
integer(i4),                                   intent(in)    :: gridcount
integer(i4),                   dimension(2),   intent(in)    :: srt
integer(i4),                   dimension(2),   intent(in)    :: cnt
integer(i4),                                   intent(in)    :: cntt
integer(i4),                                   intent(in)    :: p0
integer(i4),                                   intent(in)    :: p1
integer(i4),                                   intent(inout) :: xlen
integer(i4),                                   intent(inout) :: ylen
integer(i4),                                   intent(inout) :: tlen
real(dp),         allocatable, dimension(:),   intent(inout) :: inlon
real(dp),         allocatable, dimension(:),   intent(inout) :: inlat
real(dp),         allocatable, dimension(:),   intent(inout) :: time
integer(i4),      allocatable, dimension(:,:), intent(inout) :: indx
real(sp),                      dimension(:),   intent(inout) :: elev
type(in_metvars),              dimension(:),   intent(inout) :: invars

! Monthly input driver variables
real(sp), allocatable, dimension(:,:) :: tmp        ! mean monthly temperature (degC)
real(sp), allocatable, dimension(:,:) :: dtr        ! mean monthly diurnal temperature range (degC)
real(sp), allocatable, dimension(:,:) :: pre        ! total monthly precipitation (mm)
real(sp), allocatable, dimension(:,:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
real(sp), allocatable, dimension(:,:) :: cld        ! mean monthly cloud cover (fraction)
real(sp), allocatable, dimension(:,:) :: wnd        ! mean monthly 10m windspeed (m s-1)

integer(i4) :: ilen
integer(i4), allocatable, dimension(:) :: var_in
real(sp)    :: scale_factor
real(sp)    :: add_offset
integer(i2) :: missing_value

integer :: ifid
integer :: dimid
integer :: varid
integer :: i

! -----------------------------------------------------
! INPUT: Read dimension IDs and lengths of dimensions

ncstat = nf90_open(infile,nf90_nowrite,ifid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           ! Open netCDF-file (inpput file name, no writing rights, assigned file number)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! Check for errors (after every step)

! ncstat = nf90_inq_dimid(ifid,'lon',dimid)              ! get dimension ID from dimension 'lon' in the input file
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)   ! get dimension name and length from input file for dimension previously inquired
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inq_dimid(ifid,'lat',dimid)              ! get dimension ID from dimension 'lon' in the input file
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)   ! get dimension name and length from input file for dimension previously inquired
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
ncstat = nf90_inq_dimid(ifid,'index',dimid)             ! get dimension ID from dimension 'index' in the input file
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ilen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inq_dimid(ifid,'time',dimid)             ! Get dimension ID for time
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)   ! Get length of dimension 'time' and assign it to variable tlen
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ! ----------------------------------------------------
! ! Read variable IDs and values
!
! allocate(inlon(xlen))         ! Allocate length to longitude array
! allocate(inlat(ylen))         ! Allocate length to latitude array
! allocate(time(tlen))        ! Allocate length to time array
! allocate(indx(xlen,ylen))   ! Allocate length to index array
!
! ncstat = nf90_inq_varid(ifid,"lon",varid)                ! Get variable ID for longitude
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(ifid,varid,inlon)                  ! Get variable values for longitude
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inq_varid(ifid,"lat",varid)                ! Get variable ID for latitude
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(ifid,varid,inlat)                  ! Get variable values for latitude
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inq_varid(ifid,"index",varid)              ! Get variable ID for index
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(ifid,varid,indx)                   ! Get variable values for index
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inq_varid(ifid,"time",varid)               ! Get variable ID for time
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(ifid,varid,time)                   ! Get variable values for time
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----------------------------------------------------
! Read in elevation data from input file

allocate(var_in(ilen))

ncstat = nf90_inq_varid(ifid,'elv',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value)
  elev = real(var_in) * scale_factor + add_offset
! elsewhere
!   elev = -9999.
end where

! ----------------------------------------------------
! Read in variables from netcdf file by each gridcell
allocate(tmp(gridcount,cntt))
allocate(dtr(gridcount,cntt))
allocate(pre(gridcount,cntt))
allocate(wet(gridcount,cntt))
allocate(cld(gridcount,cntt))
allocate(wnd(gridcount,cntt))

call readdata(ifid,'tmp',srt,cnt,tmp(:,p0:p1))
call readdata(ifid,'dtr',srt,cnt,dtr(:,p0:p1))
call readdata(ifid,'pre',srt,cnt,pre(:,p0:p1))
call readdata(ifid,'wet',srt,cnt,wet(:,p0:p1))
call readdata(ifid,'cld',srt,cnt,cld(:,p0:p1))
call readdata(ifid,'wnd',srt,cnt,wnd(:,p0:p1))

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! write(0,*) 'ok read data', xlen, ylen, tlen, ilen, cnt

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

  invars(i)%tmp(:) = tmp(i,:)
  invars(i)%dtr(:) = dtr(i,:)
  invars(i)%pre(:) = pre(i,:)
  invars(i)%wet(:) = wet(i,:)
  invars(i)%cld(:) = cld(i,:)
  invars(i)%wnd(:) = wnd(i,:)

end do

end subroutine metdatainput

!---------------------------------------------------------------------

subroutine soildatainput(soilfile,gridstart,gridcount,soilvars)

! Read in the full array all gridcells and all layers, of soil variales (sand, clay and cfvo composition)

use parametersmod, only : i2,i4,sp,dp
use statevarsmod,  only : soil_vars
use statevarsmod,  only : nl

implicit none

character(100),                intent(in)    :: soilfile
integer(i4),                   intent(in)    :: gridstart
integer(i4),                   intent(in)    :: gridcount
type(soil_vars), dimension(:), intent(inout) :: soilvars

! Local variables
integer(i2), allocatable, dimension(:,:) :: var_in
real(sp),    allocatable, dimension(:,:) :: values

real(sp)    :: scale_factor
real(sp)    :: add_offset
integer(i2) :: missing_value

integer :: sfid
integer :: dimid
integer :: varid
integer :: i

! -----------------------------------------------------
! INPUT: Read in layer depth variable (zpos)

ncstat = nf90_open(soilfile,nf90_nowrite,sfid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           !
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_inq_varid(sfid,"zpos",varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(sfid,varid,zpos)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

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
! INPUT: Read in clay composition

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
! INPUT: Read in course gravel composition

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

!---------------------------------------------------------------------

subroutine topodatainput(topofile,gridstart,gridcount,topovars,validcell)

! Read in the full array all gridcells and all layers, of topographic variales (slope and aspect)

use statevarsmod, only : topo_vars

implicit none

character(100),                intent(in)    :: topofile
integer(i4),                   intent(in)    :: gridstart
integer(i4),                   intent(in)    :: gridcount
type(topo_vars), dimension(:), intent(inout) :: topovars
logical,         dimension(:), intent(inout) :: validcell

integer(i2), allocatable, dimension(:) :: var_in
real(sp),    allocatable, dimension(:) :: var_in_r
real(sp),    allocatable, dimension(:) :: values

real(sp)    :: scale_factor
real(sp)    :: add_offsetworld_slope_1800arcsec
integer(i2) :: missing_value

integer :: tfid
integer :: dimid
integer :: varid
integer :: i

! -----------------------------------------------------
! OPEN file

ncstat = nf90_open(topofile,nf90_nowrite,tfid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           !
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! -----------------------------------------------------
! INPUT: Read in elevation

allocate(var_in_r(gridcount))
allocate(values(gridcount))

ncstat = nf90_inq_varid(tfid,"elev",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(tfid,varid,var_in_r,start=[gridstart],count=[gridcount])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(tfid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_get_att(tfid,varid,"scale_factor",scale_factor)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in_r /= missing_value)
  values = var_in_r !* scale_factor
elsewhere
  values    = -9999.
  validcell = .false.
end where

do i = 1, gridcount
  topovars(i)%elev = values(i)
end do

deallocate(var_in_r)
deallocate(values)

! -----------------------------------------------------
! INPUT: Read in slope

allocate(var_in_r(gridcount))
allocate(values(gridcount))

ncstat = nf90_inq_varid(tfid,"slope",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(tfid,varid,var_in_r,start=[gridstart],count=[gridcount])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(tfid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in_r /= missing_value)
  values = var_in_r
elsewhere
  values = -9999.
end where

do i = 1, gridcount
  topovars(i)%slope = values(i)
end do

deallocate(var_in_r)
deallocate(values)

! -----------------------------------------------------
! INPUT: Read in aspect

allocate(var_in_r(gridcount))
allocate(values(gridcount))

ncstat = nf90_inq_varid(tfid,"aspect",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(tfid,varid,var_in_r,start=[gridstart],count=[gridcount])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(tfid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in_r /= missing_value)
  values = var_in_r
elsewhere
  values = -9999.
end where

do i = 1, gridcount
  topovars(i)%aspect = values(i)
end do

deallocate(var_in_r)
deallocate(values)

! -----------------------------------------------------
! INPUT: Read in gridcell area

allocate(var_in_r(gridcount))
allocate(values(gridcount))

ncstat = nf90_inq_varid(tfid,"cell_area",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(tfid,varid,var_in_r,start=[gridstart],count=[gridcount])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(tfid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in_r /= missing_value)
  values = var_in_r
elsewhere
  values = -9999.
end where

do i = 1, gridcount
  topovars(i)%cellarea = values(i)
end do

deallocate(var_in_r)
deallocate(values)

! -----------------------------------------------------
! INPUT: Read in areafrac

allocate(var_in_r(gridcount))
allocate(values(gridcount))

ncstat = nf90_inq_varid(tfid,"areafrac",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(tfid,varid,var_in_r,start=[gridstart],count=[gridcount])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(tfid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in_r /= missing_value)
  values = var_in_r
elsewhere
  values = -9999.
end where

do i = 1, gridcount
  topovars(i)%areafrac = values(i)
end do

deallocate(var_in_r)
deallocate(values)

!------

ncstat = nf90_close(tfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


end subroutine topodatainput

!---------------------------------------------------------------------

subroutine LAIdatainput(LAIfile,gridcount,srt,cnt,ndyear,vegvars)

! Read in the full array all gridcells and all layers, of soil variales (sand, clay and cfvo composition)

use statevarsmod, only : veg_vars
use utilitiesmod, only : getmonth

implicit none

character(100),                intent(in)    :: LAIfile
integer(i4),                   intent(in)    :: gridcount
integer(i4),     dimension(2), intent(in)    :: srt
integer(i4),     dimension(2), intent(in)    :: cnt
integer(i4),                   intent(in)    :: ndyear
type(veg_vars),  dimension(:), intent(inout) :: vegvars

integer(i4) :: gridstart

real(sp), allocatable, dimension(:,:) :: var_in
real(sp), allocatable, dimension(:,:) :: values

real(sp)    :: scale_factor
real(sp)    :: add_offset
integer(i2) :: missing_value

integer(i4) :: month
integer(i4) :: sday
integer(i4) :: eday

integer :: ifid
integer :: dimid
integer :: varid
integer :: i
integer :: d

!------

gridstart = srt(1)

! -----------------------------------------------------
! INPUT: Read in layer depth variable (zpos)

ncstat = nf90_open(LAIfile,nf90_nowrite,ifid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           !
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


! -----------------------------------------------------
! INPUT: Read in LAI data file

allocate(var_in(gridcount,12))
allocate(values(gridcount,12))

ncstat = nf90_inq_varid(ifid,"LAI",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[gridstart,1],count=[gridcount,12])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

scale_factor = 1.0

where (var_in /= missing_value)
  values = var_in
elsewhere
  values = -9999.
end where

do i = 1, gridcount
  do d = 1, 366

    call getmonth(d,ndyear,month,sday,eday)

    ! vegvars(i,d)%lai = values(i,month)

  end do
end do

deallocate(var_in)
deallocate(values)


!------

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


end subroutine LAIdatainput

!---------------------------------------------------------------------

end module netcdfinputmod
