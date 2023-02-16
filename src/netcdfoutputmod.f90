module netcdfoutputmod

! use outputmod,     only : mpivars
! use pftparmod,     only : npft,pftpar,tree
! use statevarsmod,  only : dayvars,soilvars,gppvars,vegvars,ndyear,nd,calcyrs,ilen
! use biome1mod,     only : biomevars,biome_mean
! use weathergenmod, only : roundto


implicit none

!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

subroutine netcdfoutput(outfile,ssrt,ccnt,year,ndyear,nd,sv)

use parametersmod, only : i2,i4,sp
use errormod,      only : ncstat,netcdf_err
use mpivarsmod,    only : ofid
use pftparmod,     only : npft,tree
use statevarsmod,  only : statevars
use netcdf
use mpi

implicit none

character(100),                intent(in) :: outfile
integer(i4),     dimension(2), intent(in) :: ssrt
integer(i4),     dimension(2), intent(in) :: ccnt
integer(i4),                   intent(in) :: year
integer(i4),                   intent(in) :: ndyear
integer(i4),     dimension(:), intent(in) :: nd
type(statevars), dimension(:), intent(in) :: sv

! integer(i4) :: ofid
integer(i4) :: dimid
integer(i4) :: varid
integer(i4), dimension(2) :: srt
integer(i4), dimension(2) :: cnt

integer(i2), allocatable, dimension(:) :: outvar_biome
integer(i2), allocatable, dimension(:) :: outvar
real(sp),    allocatable, dimension(:) :: outvar_r
real(sp),    allocatable, dimension(:,:,:) :: outvar_pft

real(sp), dimension(npft) :: fpc_grid
real(sp), dimension(npft) :: meanfpc

integer(i4) :: i
integer(i4) :: pft
integer(i4) :: grid

!-------------------------

srt(1) = ssrt(1)
cnt(1) = ccnt(1)

if (year == 1) then
  srt(2) = 1
else
  srt(2) = sum(nd(1:(year-1)*12)) + 1
end if

cnt(2) = ndyear

!-------------------------

allocate(outvar_biome(cnt(1)))

allocate(outvar(cnt(1)))

allocate(outvar_r(cnt(1)))

allocate(outvar_pft(cnt(1),1,npft))

!-------------------------

if (year == 1) then
  ncstat = nf90_open(outfile,nf90_write,ofid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
end if

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
!     outvar_r(grid,i) = sum(sv(grid)%gppvars%gpp)
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
!     outvar_r(grid,i) = sum(sv(grid)%gppvars%npp)
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

ncstat = nf90_inq_varid(ofid,'gpp',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do grid = 1, cnt(1)
  do pft = 1, npft
    outvar_pft(grid,1,pft) = sum(sv(grid)%gppvars%gpp(:,pft))
  end do
end do

ncstat = nf90_put_var(ofid,varid,outvar_pft,start=[srt(1),year,1],count=[cnt(1),1,npft])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

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

! ncstat = nf90_close(ofid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


end subroutine netcdfoutput

!---------------------------------------------------------------------

subroutine netcdfoutput_onelayer(outfile,ssrt,ccnt,year,ndyear,nd,sv)

use parametersmod, only : i2,i4,sp
use errormod,      only : ncstat,netcdf_err
use pftparmod,     only : npft,tree
use statevarsmod,  only : statevars
use netcdf
use mpi

implicit none

character(100),                intent(in) :: outfile
integer(i4),     dimension(2), intent(in) :: ssrt
integer(i4),     dimension(2), intent(in) :: ccnt
integer(i4),                   intent(in) :: year
integer(i4),                   intent(in) :: ndyear
integer(i4),     dimension(:), intent(in) :: nd
type(statevars), dimension(:), intent(in) :: sv

integer(i4) :: ofid
integer(i4) :: dimid
integer(i4) :: varid
integer(i4), dimension(2) :: srt
integer(i4), dimension(2) :: cnt

integer(i2), allocatable, dimension(:)   :: outvar
real(sp),    allocatable, dimension(:)   :: outvar_r
real(sp),    allocatable, dimension(:,:) :: outvar_pft

real(sp), dimension(npft) :: fpc_grid
real(sp), dimension(npft) :: meanfpc

integer(i4) :: i
integer(i4) :: grid

!-------------------------

srt(1) = ssrt(1)
cnt(1) = ccnt(1)

if (year == 1) then
  srt(2) = 1
else
  srt(2) = sum(nd(1:(year-1)*12)) + 1
end if

cnt(2) = ndyear

!-------------------------

allocate(outvar(cnt(1)))

allocate(outvar_r(cnt(1)))

allocate(outvar_pft(cnt(1),npft))

!-------------------------

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
    outvar_r(grid) = sum(sv(grid)%dayvars%rhum(1:ndyear)) / ndyear
    if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'dpet',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
      outvar_r(grid) = sum(sv(grid)%dayvars%dpet(1:ndyear))
      if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'daet',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
      outvar_r(grid) = sum(sv(grid)%dayvars%daet(1:ndyear))
      if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'alpha',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    outvar_r(grid) = sum(sv(grid)%dayvars%alpha(1:ndyear)) / ndyear
    if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'gpp',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    ! outvar_r(grid) = sum(sv(grid)%gppvars%gpp(:,1)) + sum(sv(grid)%gppvars%gpp(:,2)) + sum(sv(grid)%gppvars%gpp(:,6)) &
    !                   + sum(sv(grid)%gppvars%gpp(:,3)) + sum(sv(grid)%gppvars%gpp(:,4)) + sum(sv(grid)%gppvars%gpp(:,5)) &
    !                   + sum(sv(grid)%gppvars%gpp(:,7)) + sum(sv(grid)%gppvars%gpp(:,8))

    outvar_r(grid) = sv(grid)%longave%gpp

    if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'npp',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    ! outvar_r(grid) = sum(sv(grid)%gppvars%npp(:,1)) + sum(sv(grid)%gppvars%npp(:,2)) + sum(sv(grid)%gppvars%npp(:,6)) &
    !                   + sum(sv(grid)%gppvars%npp(:,3)) + sum(sv(grid)%gppvars%npp(:,4)) + sum(sv(grid)%gppvars%npp(:,5)) &
    !                   + sum(sv(grid)%gppvars%npp(:,7)) + sum(sv(grid)%gppvars%npp(:,8))

    outvar_r(grid) = sv(grid)%longave%npp

    if (sv(grid)%gppvars%npp(1,1) == -9999.) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'lm_ind',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    ! outvar_r(grid) = sum(sv(grid)%vegvars%lm_ind(:))
    ! outvar_r(grid) = sum(sv(grid)%soilvars%Tliq(1:6)) / sum(sv(grid)%soilvars%Tsat(1:6))
    outvar_r(grid) = sv(grid)%longave%soilw
    if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'rm_ind',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    ! outvar_r(grid) = sum(sv(grid)%vegvars%rm_ind(:))
    ! outvar_r(grid) = sv(grid)%soilvars%zsno
    ! outvar_r(grid) = sum(sv(grid)%soilvars%Tliq(1:6)) + sum(sv(grid)%soilvars%Tice(1:6)) / sum(sv(grid)%soilvars%Tsat(1:6))
    ! outvar_r(grid) = sum(sv(grid)%soilvars%Tsoil(1:3))/3. - 273.15
    ! outvar_r(grid) = sv(grid)%vegvars%rm_ind(1)
    outvar_r(grid) = sv(grid)%longave%Tsoil - 273.15
    if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'sm_ind',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    ! outvar_r(grid) = sum(sv(grid)%vegvars%lm_ind(:) * sv(grid)%vegvars%nind(:)) +&
    !                   sum(sv(grid)%vegvars%rm_ind(:) * sv(grid)%vegvars%nind(:))
    outvar_r(grid) = sv(grid)%longave%lm_ind + sv(grid)%longave%rm_ind + sv(grid)%longave%sm_ind
    if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'hm_ind',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    ! outvar_r(grid) =  sum(sv(grid)%vegvars%sm_ind(:) * sv(grid)%vegvars%nind(:)) +&
    !                   sum(sv(grid)%vegvars%hm_ind(:) * sv(grid)%vegvars%nind(:))
    ! outvar_r(grid) = sv(grid)%soilvars%zsno
    outvar_r(grid) = sv(grid)%longave%zsno
    if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'height',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)
    if (count(sv(grid)%vegvars%present .and. tree(1:npft)) /= 0) then
      ! outvar_r(grid) = sum(sv(grid)%vegvars%height) / count(sv(grid)%vegvars%present .and. tree(1:npft))
      ! outvar_r(grid) = maxval(sv(grid)%vegvars%height)
      outvar_r(grid) = sv(grid)%longave%height
    else
      outvar_r(grid) = 0.
    end if

    if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'fpc_max',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)

    fpc_grid = sv(grid)%vegvars%fpc_grid
    meanfpc  = sv(grid)%vegvars%meanfpc

    if (meanfpc(1) == maxval(meanfpc)) outvar_r(grid) = 1.
    if (meanfpc(2) == maxval(meanfpc)) outvar_r(grid) = 2.
    if (meanfpc(3) == maxval(meanfpc)) outvar_r(grid) = 3.
    if (meanfpc(4) == maxval(meanfpc)) outvar_r(grid) = 4.
    if (meanfpc(5) == maxval(meanfpc)) outvar_r(grid) = 5.
    if (meanfpc(6) == maxval(meanfpc)) outvar_r(grid) = 6.
    if (meanfpc(7) == maxval(meanfpc)) outvar_r(grid) = 7.
    if (meanfpc(8) == maxval(meanfpc)) outvar_r(grid) = 8.
    if (meanfpc(9) == maxval(meanfpc)) outvar_r(grid) = 9.

    if (sum(meanfpc) < 0.1) outvar_r(grid) = 0.

    if (.not.(sv(grid)%validcell)) outvar_r(grid) = -9999.

  end do

ncstat = nf90_put_var(ofid,varid,outvar_r,start=[srt(1)],count=[cnt(1)])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! -----------------------------------------------------------

ncstat = nf90_inq_varid(ofid,'cover',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do grid = 1, cnt(1)

    outvar_pft(grid,1:npft) = sv(grid)%vegvars%meanfpc

    if (.not.(sv(grid)%validcell)) outvar_pft(grid,:) = -9999.
  end do

ncstat = nf90_put_var(ofid,varid,outvar_pft,start=[srt(1),1],count=[cnt(1),npft])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

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
