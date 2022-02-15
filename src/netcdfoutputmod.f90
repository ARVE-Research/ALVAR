module netcdfoutputmod

use parametersmod, only : i2,i4,sp,dp,missing_i2
use errormod,      only : ncstat,netcdf_err
use outputmod,     only : infompi
use pftparmod,     only : npft,pftpar,tree
use statevarsmod,  only : dayvars,soilvars,gppvars,vegvars,ndyear,nd,calcyrs,ilen
use biome1mod,     only : biomevars,biome_mean
use weathergenmod, only : roundto
use netcdf
use mpi

implicit none

!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

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

integer :: i,grid



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
! ncstat = nf90_inq_varid(ofid,'biome1_year',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_biome = biomevars(:,year)%biome
!
! where (outvar_biome == missing_i2)
!   outvar_biome = missing_i2
! else where
!   outvar_biome = nint(outvar_biome / 0.1)
! end where
!
! ncstat = nf90_put_var(ofid,varid,outvar_biome,start=[srt(1),year],count=[cnt(1),1])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! !-----------------------------------------------------------
!
! if (year == calcyrs) then
!
!   ncstat = nf90_inq_varid(ofid,'biome1_mean',varid)
!   if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
!   outvar_biome = biome_mean
!
!   where (outvar_biome == missing_i2)
!     outvar_biome = missing_i2
!   else where
!     outvar_biome = nint(outvar_biome / 0.1)
!   end where
!
!   ncstat = nf90_put_var(ofid,varid,outvar_biome,start=[srt(1)],count=[cnt(1)])
!   if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! end if
!
! ! -----------------------------------------------------------
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
! ncstat = nf90_inq_varid(ofid,'tmean',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar = nint(dayvars(:,1:ndyear)%tmean / 0.1)
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
! -----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'wind',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar = nint(dayvars(:,1:ndyear)%wind / 0.1)
!
! ncstat = nf90_put_var(ofid,varid,outvar,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! -----------------------------------------------------------
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
!-----------------------------------------------------------

! ncstat = nf90_inq_varid(ofid,'dpet',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ! outvar_r = dayvars(:,1:ndyear)%dpet
!
! do grid = 1, cnt(1)
!   do i = 1, ndyear
!     outvar_r(grid,i) = sum(gppvars(grid,:)%gpp)
!     if (gppvars(grid,1)%gpp == -9999.) outvar_r(grid,i) = -9999.
!   end do
! end do
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! !-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'daet',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ! outvar_r = dayvars(:,1:ndyear)%daet
!
! do grid = 1, cnt(1)
!   do i = 1, ndyear
!     outvar_r(grid,i) = sum(gppvars(grid,:)%npp)
!     if (gppvars(grid,1)%npp == -9999.) outvar_r(grid,i) = -9999.
!   end do
! end do
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

! ncstat = nf90_inq_varid(ofid,'alpha',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%alpha
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! -----------------------------------------------------------

! ncstat = nf90_inq_varid(ofid,'gpp',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! do grid = 1, cnt(1)
!   do i = 1, ndyear
!     outvar_r(grid,i) = sum(vegvars(grid,:)%gpp0)
!   end do
! end do
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'ForFireMk5',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%ForFireMk5
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


end subroutine netcdfoutput

!---------------------------------------------------------------------

subroutine netcdfoutput_onelayer(info,job,year)

use pftparmod, only : npft
use statevarsmod, only : vegvars

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
integer(i2), allocatable, dimension(:) :: outvar
real(sp),    allocatable, dimension(:) :: outvar_r

real(sp), dimension(npft) :: fpc_grid

integer :: i,grid



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

allocate(outvar(cnt(1)))

allocate(outvar_r(cnt(1)))

ncstat = nf90_open(outfile,nf90_write,ofid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'biome1_year',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_biome = biomevars(:,year)%biome
!
! where (outvar_biome == missing_i2)
!   outvar_biome = missing_i2
! else where
!   outvar_biome = nint(outvar_biome / 0.1)
! end where
!
! ncstat = nf90_put_var(ofid,varid,outvar_biome,start=[srt(1),year],count=[cnt(1),1])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! !-----------------------------------------------------------
!
! if (year == calcyrs) then
!
!   ncstat = nf90_inq_varid(ofid,'biome1_mean',varid)
!   if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
!   outvar_biome = biome_mean
!
!   where (outvar_biome == missing_i2)
!     outvar_biome = missing_i2
!   else where
!     outvar_biome = nint(outvar_biome / 0.1)
!   end where
!
!   ncstat = nf90_put_var(ofid,varid,outvar_biome,start=[srt(1)],count=[cnt(1)])
!   if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! end if
!
! ! -----------------------------------------------------------
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
! ncstat = nf90_inq_varid(ofid,'tmean',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar = nint(dayvars(:,1:ndyear)%tmean / 0.1)
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
! -----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'wind',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar = nint(dayvars(:,1:ndyear)%wind / 0.1)
!
! ncstat = nf90_put_var(ofid,varid,outvar,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! -----------------------------------------------------------
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

ncstat = nf90_inq_varid(ofid,'rhum',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    outvar_r(grid) = sum(dayvars(grid,:)%rhum) / ndyear
    if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'dpet',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
      outvar_r(grid) = sum(dayvars(grid,:)%dpet)
      if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'daet',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
      outvar_r(grid) = sum(dayvars(grid,:)%daet)
      if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'alpha',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    outvar_r(grid) = sum(dayvars(grid,:)%alpha) / ndyear
    if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'gpp',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    outvar_r(grid) = sum(gppvars(grid,:)%gpp(1)) + sum(gppvars(grid,:)%gpp(2)) &
                      + sum(gppvars(grid,:)%gpp(3)) + sum(gppvars(grid,:)%gpp(4)) + sum(gppvars(grid,:)%gpp(5))
    if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'npp',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    outvar_r(grid) = sum(gppvars(grid,:)%npp(1)) + sum(gppvars(grid,:)%npp(2)) &
                      + sum(gppvars(grid,:)%npp(3)) + sum(gppvars(grid,:)%npp(4)) + sum(gppvars(grid,:)%npp(5))
    if (gppvars(grid,1)%npp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'lm_ind',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    ! outvar_r(grid) = sum(vegvars(grid)%lm_ind(:))
    outvar_r(grid) = vegvars(grid)%lm_ind(4)
    if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'rm_ind',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    ! outvar_r(grid) = sum(vegvars(grid)%rm_ind(:))
    outvar_r(grid) = vegvars(grid)%rm_ind(5)
    if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'sm_ind',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    outvar_r(grid) = sum(vegvars(grid)%sm_ind(:))
    if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'hm_ind',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    outvar_r(grid) = sum(vegvars(grid)%hm_ind(:))
    if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'height',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    if (count(vegvars(grid)%present .and. tree(1:npft)) /= 0) then
      outvar_r(grid) = sum(vegvars(grid)%height) / count(vegvars(grid)%present .and. tree(1:npft))
    else
      outvar_r(grid) = 0.
    end if

    if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'fpc_max',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)

    fpc_grid = vegvars(grid)%fpc_grid

    if (fpc_grid(1) == maxval(fpc_grid)) outvar_r(grid) = 1.
    if (fpc_grid(2) == maxval(fpc_grid)) outvar_r(grid) = 2.
    if (fpc_grid(3) == maxval(fpc_grid)) outvar_r(grid) = 3.
    if (fpc_grid(4) == maxval(fpc_grid)) outvar_r(grid) = 4.
    if (fpc_grid(5) == maxval(fpc_grid)) outvar_r(grid) = 5.

    if (sum(fpc_grid) < 0.1) outvar_r(grid) = 0.

    if (gppvars(grid,1)%gpp(1) == -9999.) outvar_r(grid) = -9999.

  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! -----------------------------------------------------------

! ncstat = nf90_inq_varid(ofid,'gpp',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! do grid = 1, cnt(1)
!   do i = 1, ndyear
!     outvar_r(grid,i) = sum(vegvars(grid,:)%gpp0)
!   end do
! end do
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------
!
! ncstat = nf90_inq_varid(ofid,'ForFireMk5',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! outvar_r = dayvars(:,1:ndyear)%ForFireMk5
!
! ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt],count=[cnt])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


end subroutine netcdfoutput_onelayer

!---------------------------------------------------------------------




end module netcdfoutputmod
