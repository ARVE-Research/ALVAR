module netcdfsetupmod

! use parametersmod, only : maxoutvars
use parametersmod, only : i4,sp,stdout,stderr

implicit none

public  :: netcdf_create
private :: getvarinfo
private :: declvar

integer(i4), parameter :: maxoutvars = 100
integer(i4), parameter :: max_ndims  = 10

type variableinfo
  character(20) :: name
  character(80) :: longname
  character(20) :: units
  integer(i4)   :: ndims
  character(20), dimension(max_ndims) :: dimnames
end type variableinfo

type(variableinfo), target, dimension(maxoutvars) :: varinfo

character(50), dimension(maxoutvars) :: reqvarsout

! character(40), dimension(maxoutvars) :: outputvar = 'null'


contains

!------------------------------------------------------------------------------------------------------------------------------

subroutine genoutfile(info)

! Subroutine to generate the list-formatted outfile name for the list run output

use mpivarsmod, only : mpivars

implicit none

type(mpivars), intent(inout) :: info

integer(i4) :: n

info%outfile = trim(info%outfile)

n = len(trim(info%outfile))

info%outfile_list = info%outfile(1:n-3)//'_list-formmatted.nc'

end subroutine genoutfile

!------------------------------------------------------------------------------------------------------------------------------

subroutine netcdf_create(info,grid)

use parametersmod, only : i4,sp,dp,netcdf_baseyr
use pftparmod,     only : npft
use errormod,      only : ncstat,netcdf_err
use mpivarsmod,    only : mpivars
use typesizes
use netcdf

! use iovariablesmod, only : ofid,lonvect,latvect,srtx,cntx,endx,srty,cnty,endy,outputfile,outputvar,cellindex,cellmask,  &
!                            calcforagers,gridres

implicit none

type(mpivars), target,   intent(inout) :: info
logical,       optional, intent(in)    :: grid

character(200) :: infile
character(200) :: outfile

integer(i4) :: xlen
integer(i4) :: ylen
integer(i4) :: ilen
integer(i4) :: tlen
integer(i4) :: ncells
integer(i4) :: tran_startyr

real(dp),    allocatable, dimension(:)   :: lon
real(dp),    allocatable, dimension(:)   :: lat
integer(i4), allocatable, dimension(:,:) :: indx
integer(i4), allocatable, dimension(:)   :: time

integer(i4), dimension(2) :: dimids

integer(i4) :: ifid
integer(i4) :: ofid
integer(i4) :: dimid
integer(i4) :: varid
integer(i4) :: i
integer(i4) :: x
integer(i4) :: y

character(40),          dimension(maxoutvars) :: outputvar
logical,      pointer, dimension(:) :: outvar_on

integer(i4) :: reqoutvars

character(8)  :: today
character(10) :: now

real(dp), dimension(2) :: xrange
real(dp), dimension(2) :: yrange

real(dp) :: xres
real(dp) :: yres

integer(i4), allocatable, dimension(:) :: pftnum

character(40), dimension(2) :: varlabel !names of variables that will be automatically output

!character(40), dimension(2), parameter :: varlabel = ['    landf' ,'foragerPD']  !names of variables that will be automatically output

!----------

infile       = info%cfile_spinup
ncells       = info%ilen
tlen         = info%transientyears
tran_startyr = info%tran_startyr

outputvar =  info%outputvar
outvar_on => info%outvar_on

!----------

outvar_on = .false.

!----------

if (grid) then
  call getvarinfo('./input/outvarsinfo.namelist')
  outfile = info%outfile

else
  call getvarinfo('./input/outvarsinfo_indx.namelist')
  outfile = info%outfile_list
end if

! xres = gridres(1)
! yres = gridres(2)
!
! write(stdout,*)'input grid resolution',xres,yres,' degrees'
!
! xrange(1) = minval(lonvect(srtx:endx)) - 0.5 * xres
! xrange(2) = maxval(lonvect(srtx:endx)) + 0.5 * xres
!
! yrange(1) = minval(latvect(srty:endy)) - 0.5 * yres
! yrange(2) = maxval(latvect(srty:endy)) + 0.5 * yres
!
! allocate(cellmask(cntx,cnty))
!
! cellmask = .false.

!----------------------------------
!dimensions and dimension variables

write(0,*)
write(0,*) 'Reading lon, lat and index from infile'

ncstat = nf90_open(infile,nf90_nowrite,ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lon',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lat',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'index',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ilen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!----------------------------------

allocate(lon(xlen))
allocate(lat(ylen))
allocate(indx(xlen,ylen))

ncstat = nf90_inq_varid(ifid,"lon",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lon)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lat)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"index",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,indx)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

xres = lon(2) - lon(1)
yres = lat(2) - lat(1)

xrange(1) = minval(lon) - 0.5 * xres
xrange(2) = maxval(lon) + 0.5 * xres

yrange(1) = minval(lat) - 0.5 * yres
yrange(2) = maxval(lat) + 0.5 * yres


write(0,*) 'Got lon, lat and index from infile'
write(0,*)

!----------------------------------

write(stdout,'(a,a)')'output file: ',trim(outfile)

ncstat = nf90_create(outfile,nf90_hdf5,ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'title','ALVAR netCDF output file')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

call date_and_time(today,now)

ncstat = nf90_put_att(ofid,nf90_global,'timestamp',today//' '//now(1:4))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'Conventions','COARDS')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'node_offset',1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

write(stdout,*)'added global atts'

!----

ncstat = nf90_def_dim(ofid,'lon',xlen,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(1) = dimid

ncstat = nf90_def_var(ofid,'lon',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','longitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_east')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',xrange)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'lat',ylen,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(2) = dimid

ncstat = nf90_def_var(ofid,'lat',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','latitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',yrange)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! index

ncstat = nf90_def_dim(ofid,'index',ilen,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'index',nf90_int,dimids(1:2),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,indx)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','index of lon and lat')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','1 to validcells')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',-1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'layer',6,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'layer',nf90_short,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','soil layer')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','layer')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

! ncstat = nf90_def_dim(ofid,'tile',size(lutype),dimid)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_def_var(ofid,'tile',nf90_short,dimid,varid)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_put_att(ofid,varid,'long_name','land use tile')
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_put_att(ofid,varid,'units','layer')
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'pft',npft,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'pft',nf90_short,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','Plant Functional Type')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','PFT')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'month',12,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'month',nf90_short,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','month')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','month')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'time',tlen,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'time',nf90_int,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','time')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','years since 1950-00-00 00:00')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!---------------------------------------------------------------------------
!regular variables

varlabel(1) =  'landf'
varlabel(2) =  'foragerPD'
outputvar = eoshift(outputvar,-1,varlabel(1))

! if (calcforagers) then
!   outputvar = eoshift(outputvar,-1,varlabel(2))
! end if

reqoutvars = count(outputvar /= 'null')

if (reqoutvars < 1) then
  write(stdout,*)'WARNING, no variables were specified for output!'

else
  write(stdout,*)'the following variables will be written to netCDF output'

  do i = 1,reqoutvars

    ! write(stdout,*) outputvar(i)
    call declvar(ofid,outputvar(i),outvar_on)

  end do
end if
!---------------------------------------------------------------------------

ncstat = nf90_enddef(ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

!write the dimension variables (except for time)

allocate(pftnum(npft))
allocate(time(tlen))

forall (i=1:npft)
  pftnum(i) = i
end forall

do i = 1, tlen
  time(i) = tran_startyr-netcdf_baseyr + (i-1)
end do

ncstat = nf90_inq_varid(ofid,'lon',varid)
ncstat = nf90_put_var(ofid,varid,lon)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'lat',varid)
ncstat = nf90_put_var(ofid,varid,lat)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'index',varid)
ncstat = nf90_put_var(ofid,varid,indx)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'layer',varid)
ncstat = nf90_put_var(ofid,varid,[1,2,3,4,5,6])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'pft',varid)
ncstat = nf90_put_var(ofid,varid,[1,2,3,4,5,6,7,8,9])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_inq_varid(ofid,'tile',varid)
! ncstat = nf90_put_var(ofid,varid,[1,2,3])
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'time',varid)
ncstat = nf90_put_var(ofid,varid,time)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'month',varid)
ncstat = nf90_put_var(ofid,varid,[1,2,3,4,5,6,7,8,9,10,11,12])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_close(ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!---------------------------------------------------------------------------

deallocate(pftnum)

end subroutine netcdf_create

!------------------------------------------------------------------------------------------------------------------------------

subroutine getvarinfo(varnamelist)

implicit none

character(*), intent(in) :: varnamelist

namelist / varinfolist / varinfo

open(87, file=varnamelist, status='old')

read(87,varinfolist)

close(87)

end subroutine getvarinfo

!--------------------------------------------

subroutine declvar(fid,varname,outvar_on)

use errormod, only : ncstat,netcdf_err
! use typesizes
use netcdf

implicit none

integer,                    intent(in)    :: fid
character(*),               intent(in)    :: varname
logical,      dimension(:), intent(inout) :: outvar_on

character(20), pointer :: name
character(80), pointer :: longname
character(20), pointer :: units

character(20), pointer, dimension(:) :: dimnames

real :: missing = -9999.

integer, pointer :: ndims

integer :: i,j

integer :: dimlen
integer :: varid

integer, allocatable, dimension(:) :: dimids
integer, allocatable, dimension(:) :: chunks

!------
!scan varinfo for the variable to be processed

do i = 1,maxoutvars
  if (trim(varinfo(i)%name) == trim(varname)) then
    j = i
    outvar_on(i) = .true.
    exit
  else
    j = 0
  end if
end do

if (j == 0) then
  write(stdout,*)'ERROR, a variable was requested for output that was not present'
  write(stdout,*)'in the attribute table: outvarsinfo.namelist'
  write(stdout,*)'Please check the variable name or add metadata to the table.'
  stop
end if

name     => varinfo(i)%name
longname => varinfo(i)%longname
units    => varinfo(i)%units
ndims    => varinfo(i)%ndims
dimnames => varinfo(i)%dimnames

allocate(dimids(ndims))
allocate(chunks(ndims))

do i = 1,ndims

  ncstat = nf90_inq_dimid(fid,dimnames(i),dimids(i))
  if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inquire_dimension(fid,dimids(i),len=dimlen)
  if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

  if (i <= 2) then
    chunks(i) = dimlen
  else
    chunks(i) = 1
  end if

end do

ncstat = nf90_def_var(fid,varname,nf90_float,dimids,varid) !,chunksizes=chunks,deflate_level=1,shuffle=.false.)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(fid,varid,'long_name',longname)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(fid,varid,'units',units)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(fid,varid,'_FillValue',missing)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(fid,varid,'missing_value',missing)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

deallocate(dimids)
deallocate(chunks)

end subroutine declvar

!------------------------------------------------------------------------------------------------------------------------------

end module netcdfsetupmod
