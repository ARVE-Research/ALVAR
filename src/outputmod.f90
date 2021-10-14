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

integer, dimension(5) :: dimids

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
! days

ncstat = nf90_def_dim(ofid,'days',ndays,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(4) = dimid

ncstat = nf90_def_var(ofid,'days',nf90_int,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,(/(i,i=1,ndays,1)/))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','time in number of days')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','days')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! years

ncstat = nf90_def_dim(ofid,'years',calcyrs,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(5) = dimid

ncstat = nf90_def_var(ofid,'years',nf90_int,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,(/(i,i=1,calcyrs,1)/))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','time in number of years')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','years')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)


!----
!biome type

ncstat = nf90_def_var(ofid,'biome1_year',nf90_short,[dimids(3),dimids(5)],varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','BIOME1 annual output')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','type')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'scale_factor',0.1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!biome type

ncstat = nf90_def_var(ofid,'biome1_mean',nf90_short,dimids(3),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','BIOME1 mean climate output')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','type')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'scale_factor',0.1)
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
!average daytime temperature

ncstat = nf90_def_var(ofid,'tday',nf90_short,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily average daytimetemperature')
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
!average nighttime temperature

ncstat = nf90_def_var(ofid,'tnight',nf90_short,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily average nighttime temperature')
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
!average 24h temperature

ncstat = nf90_def_var(ofid,'tmean',nf90_short,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily average 24h temperature')
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
!mean dewpoint temperature

ncstat = nf90_def_var(ofid,'tdew',nf90_short,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily mean dewpoint temperature')
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
!mean 10m windspeed

ncstat = nf90_def_var(ofid,'wind',nf90_short,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily mean 10m windspeed')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','m s-1')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'scale_factor',0.1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!24 hour total precipitation

ncstat = nf90_def_var(ofid,'prec',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily precipitation')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','mm d-1')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!daily mean relative humidity

ncstat = nf90_def_var(ofid,'dayl',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daylength')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','h')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!Pjj

ncstat = nf90_def_var(ofid,'vpd',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','vapor pressure deficit')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','Pa')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!daily mean relative humidity

ncstat = nf90_def_var(ofid,'srad',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','downwelling surface shortwave radiation')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','kJ m-2 d-1')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!daily mean relative humidity

ncstat = nf90_def_var(ofid,'rhum',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily mean relative humidity')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','%')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!daily PET

ncstat = nf90_def_var(ofid,'dpet',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily potential evapotranspiration')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','mm d-1')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!daily PET

ncstat = nf90_def_var(ofid,'daet',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','daily actual evapotranspiration')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','mm d-1')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!daily PET

ncstat = nf90_def_var(ofid,'alpha',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','AET/PET ratio')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','fraction')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!GPP

ncstat = nf90_def_var(ofid,'gpp',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','gross primary productivity')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','g C m-2 month-1')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!forest fire danger meter Mk5

ncstat = nf90_def_var(ofid,'ForFireMk5',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','Forest fire danger meter Mark 5')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','unit')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!forest fire danger meter Mk5

ncstat = nf90_def_var(ofid,'sloperad',nf90_double,dimids(3:4),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','slope to flat direct beam radiation ratio')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','ratio')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)




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
            'prec ', &
            'dayl ', &
            'dsol ', &
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
          dayvars(grid,day)%prec, &
          dayvars(grid,day)%dayl, &
          dayvars(grid,day)%dsol, &
          dayvars(grid,day)%rhum, &
          dayvars(grid,day)%srad, &
          dayvars(grid,day)%dpet, &
          dayvars(grid,day)%vpd



end subroutine printgrid


!-------------------------------------------------------------------------------------------------

end module outputmod
