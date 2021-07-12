module netcdfoutputmod

use parametersmod, only : i2,i4,sp,dp,missing_i2
use errormod,      only : ncstat,netcdf_err
use outputmod,     only : infompi
use metvarsmod,    only : dayvars,ndyear,nd,calcyrs
use simplesoilmod, only : soilvars
use biome1mod,     only : biomevars,biome_mean
use weathergenmod, only : roundto
use netcdf
use mpi

implicit none

contains

subroutine netcdfoutput(info,job,year)

type(infompi), target    , intent(in) :: info
integer(i4), dimension(2), intent(in) :: job
integer                  , intent(in) :: year

character(100), pointer :: outfile
integer :: ofid
integer :: dimid
integer :: varid

integer(i4), dimension(2) :: srt
integer(i4), dimension(2) :: cnt

integer(i2), allocatable, dimension(:)   :: outvar_biome
integer(i2), allocatable, dimension(:,:) :: outvar
real(sp),    allocatable, dimension(:,:) :: outvar_r

integer :: i



outfile => info%outfile


srt(1) = job(1)
cnt(1) = job(2)


if (year == 1) then
  srt(2) = 1
else
  srt(2) = sum(nd(1:(year-1)*12)) + 1
end if

cnt(2) = ndyear

allocate(outvar_biome(cnt(1)))

allocate(outvar(cnt(1),ndyear))

allocate(outvar_r(cnt(1),ndyear))

ncstat = nf90_open(outfile,nf90_write,ofid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------
!
ncstat = nf90_inq_varid(ofid,'biome1_year',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

outvar_biome = biomevars(:,year)%biome

where (outvar_biome == missing_i2)
  outvar_biome = missing_i2
else where
  outvar_biome = nint(outvar_biome / 0.1)
end where

ncstat = nf90_put_var(ofid,varid,outvar_biome,start=[srt(1),year],count=[cnt(1),1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

if (year == calcyrs) then

  ncstat = nf90_inq_varid(ofid,'biome1_mean',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  outvar_biome = biome_mean

  where (outvar_biome == missing_i2)
    outvar_biome = missing_i2
  else where
    outvar_biome = nint(outvar_biome / 0.1)
  end where

  ncstat = nf90_put_var(ofid,varid,outvar_biome,start=[srt(1)],count=[cnt(1)])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end if
!
!-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'tmin',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar = nint(dayvars(:,1:ndyear)%tmin / 0.1)
!
! ncstat = nf90_put_var(ofid,varid,outvar,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
!
! !-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'tmax',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar = nint(dayvars(:,1:ndyear)%tmax / 0.1)
!
! ncstat = nf90_put_var(ofid,varid,outvar,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ! -----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'tday',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar = nint(dayvars(:,1:ndyear)%tday / 0.1)
!
! ncstat = nf90_put_var(ofid,varid,outvar,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
!
! !-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'tnight',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar = nint(dayvars(:,1:ndyear)%tnight / 0.1)
!
! ncstat = nf90_put_var(ofid,varid,outvar,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! !-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'tdew',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar = nint(dayvars(:,1:ndyear)%tdew / 0.1)
!
! ncstat = nf90_put_var(ofid,varid,outvar,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ! -----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'prec',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%prec
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ! -----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'dayl',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%dayl
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ! -----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'vpd',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%vpd
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! !-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'srad',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%srad
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! !-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'rhum',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%rhum
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! !-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'dpet',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%dpet
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! !-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'daet',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%daet
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! !-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'alpha',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%alpha
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)




ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


end subroutine netcdfoutput




end module netcdfoutputmod
