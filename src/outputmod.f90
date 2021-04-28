module outputmod

use parametersmod, only : i2,i4,sp,dp,ndaymonth
use errormod,      only : ncstat,netcdf_err
use netcdf
use mpi

implicit none

public :: getoutfile
public :: putlonlat

real(sp),    parameter :: missing_sp = -9999.
integer(i2), parameter :: missing_i2 = -32768

!--------------------
! Derived type containing job information to be sent to all processes

type infompi
  character(100) :: infile          ! infile name
  character(100) :: outfile         ! outfile name
  character(100) :: timestring      ! timestring with startyr/calcyr
  integer(i4)    :: nproc           ! number of processors
  integer(i4)    :: validcell       ! number of validcells (no longer in use)
  integer(i4)    :: t0              ! start of time dimension for reading input file
  integer(i4)    :: nt              ! count of time dimension for reading input file
  real(sp)       :: plon            ! lon of printing gridcell
  real(sp)       :: plat            ! lat of printing gridcell
end type

!--------------------

contains

!-------------------------------------------------------------------------------------------------

subroutine getoutfile(info)

use metvarsmod, only : startyr,calcyrs

implicit none

type(infompi), target, intent(in) :: info

character(100), pointer :: infile
character(100), pointer :: outfile
integer(i4) ,   pointer :: validcell

integer :: ifid
integer :: ofid
integer :: ncstat
integer :: varid
integer :: dimid

integer, dimension(4) :: dimids

integer :: xlen
integer :: ylen
integer :: ilen

real(dp),    allocatable, dimension(:)   :: lon
real(dp),    allocatable, dimension(:)   :: lat
integer(i4), allocatable, dimension(:,:) :: indx

integer(i4) :: nmonths
integer(i4) :: ndays
integer(i4), allocatable, dimension(:) :: nd

integer :: yr
integer :: mon
integer :: i

character(8) :: today
character(10) :: now

!--------------------

infile    => info%infile
outfile   => info%outfile
validcell => info%validcell

!--------------------

nmonths = 12 * calcyrs

allocate(nd(nmonths))

i = 1
do yr = startyr, (startyr+calcyrs-1)
  do mon = 1, 12

    nd(i) = ndaymonth(yr,mon)

    i = i + 1

  end do
end do

ndays = sum(nd)

!--------------------

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

!------

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

write(0,*) 'Got lon, lat and index from infile'
write(0,*)

!--------------------

write(0,*) 'Creating outfile'

ncstat = nf90_create(outfile,nf90_netcdf4,ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_create_par(outfile,ior(nf90_netcdf4,nf90_mpiio),MPI_COMM_WORLD,MPI_INFO_NULL,ofid)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'title','weathergen parallel output file')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

call date_and_time(today,now)

ncstat = nf90_put_att(ofid,nf90_global,'created',today//' '//now(1:4))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'Conventions','COARDS')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'node_offset',1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! lon

ncstat = nf90_def_dim(ofid,'lon',xlen,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(1) = dimid

ncstat = nf90_def_var(ofid,'lon',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lon)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','longitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_east')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! lat

ncstat = nf90_def_dim(ofid,'lat',ylen,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(2) = dimid

ncstat = nf90_def_var(ofid,'lat',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lat)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','latitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! index

ncstat = nf90_def_dim(ofid,'index',validcell,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(3) = dimid

ncstat = nf90_def_var(ofid,'index',nf90_int,dimids(1:2),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,indx)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','index of lon and lat')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','1 to length')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',-1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! time

ncstat = nf90_def_dim(ofid,'time',ndays,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(4) = dimid

ncstat = nf90_def_var(ofid,'time',nf90_int,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,(/(i,i=1,ndays,1)/))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','time in number of dats')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','days')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)


!----
!minimum temperature

ncstat = nf90_def_var(ofid,'tmin',nf90_short,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily minimum temperature')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degC')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'scale_factor',0.1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!maximum temperature

ncstat = nf90_def_var(ofid,'tmax',nf90_short,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily maximum temperature')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degC')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'scale_factor',0.1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_enddef(ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

write(0,*) 'Finished create outfile'
write(0,*)


end subroutine getoutfile

!-------------------------------------------------------------------------------------------------

subroutine putlonlat(ofid,id,lon,lat)

use parametersmod, only : i4,sp,dp
use netcdf
use coordsmod,     only : index
use errormod,      only : ncstat,netcdf_err

implicit none

integer,                intent(in) :: ofid
type(index),            intent(in) :: id
real(dp), dimension(:), intent(in) :: lon
real(dp), dimension(:), intent(in) :: lat

integer :: varid

!---

ncstat = nf90_inq_varid(ofid,'lon',varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lon)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'lat',varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lat)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

end subroutine putlonlat

!-------------------------------------------------------------------------------------------------

subroutine printgrid(info,grid,year,day)

! Subroutine to print user-specified lon/lat (grid) from command line (Leo Lai Apr 2021)

use metvarsmod, only : startyr,dayvars

implicit none

type(infompi), intent(in) :: info
integer      , intent(in) :: grid
integer      , intent(in) :: year
integer      , intent(in) :: day


if (year == 1 .AND. day == 1) then

  print *,  'lon ', &
            'lat ', &
            'year ', &
            'day ', &
            'tmin ', &
            'tmax ', &
            'tmean ', &
            'tday ', &
            'tnight ', &
            'tdew ', &
            'rhum ', &
            'srad ', &
            'dpet ', &
            'vpd '

end if


! Print variables onto terminal
print *,  info%plon, &
          info%plat, &
          startyr+year-1, &
          day, &
          dayvars(grid,day)%tmin, &
          dayvars(grid,day)%tmax, &
          dayvars(grid,day)%tmean, &
          dayvars(grid,day)%tday, &
          dayvars(grid,day)%tnight, &
          dayvars(grid,day)%tdew, &
          dayvars(grid,day)%rhum, &
          dayvars(grid,day)%srad, &
          dayvars(grid,day)%dpet, &
          dayvars(grid,day)%vpd



end subroutine printgrid


!-------------------------------------------------------------------------------------------------

end module outputmod