module netcdfoutputmod

use parametersmod, only : sp,i2

implicit none

public  :: netcdf_output
public  :: netcdf_close
private :: putvar2d
! private :: putvar3d

real(sp),    parameter :: rmissing = -9999.
integer(i2), parameter :: imissing = -32768

contains

!----------------------------------------------------------------------------------------------------------------------------

subroutine netcdf_output(ofid,gridstart,gridcount,year,outvar_on,sv)

use parametersmod, only : i4,sp
use pftparmod, only : npft
use statevarsmod, only : nl,statevars
use errormod, only : ncstat,netcdf_err
use netcdf
use mpi

implicit none

integer(i4),                   intent(in) :: ofid
integer(i4),                   intent(in) :: gridstart
integer(i4),                   intent(in) :: gridcount
integer(i4),                   intent(in) :: year
logical,         dimension(:), intent(in) :: outvar_on
type(statevars), dimension(:), intent(in) :: sv

integer(i4), dimension(1) :: tval
integer(i4) :: tpos

integer :: varid
integer :: i
integer :: j

real(sp), allocatable, dimension(:)     :: rvar1d
real(sp), allocatable, dimension(:,:)   :: rvar2d
real(sp), allocatable, dimension(:,:,:) :: rvar3d

real(sp), dimension(3) :: NBP

!----------------------------
!write the time variable


tval = year
tpos = year

! ncstat = nf90_inq_varid(ofid,'time',varid)
! ncstat = nf90_put_var(ofid,varid,tval,start=[tpos],count=[1])
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)


!-----------------------------------------------------
!write other variables
!-----------------------------------------------------

allocate(rvar1d(gridcount))

!----------------------------
! livebiomass
if (outvar_on(2)) then

  do i = 1, gridcount
    rvar1d(i) = sv(i)%outvars%livebiomass
    if (.not.sv(i)%validcell) rvar1d(i) = rmissing
  end do

  call putvar1d(ofid,tpos,'livebiomass',gridstart,gridcount,rvar1d)

end if

!----------------------------
! GPP
if (outvar_on(13)) then

  do i = 1, gridcount
    rvar1d(i) = sv(i)%outvars%GPP
    if (.not.sv(i)%validcell) rvar1d(i) = rmissing
  end do

  call putvar1d(ofid,tpos,'GPP',gridstart,gridcount,rvar1d)

end if

!----------------------------
! NPP
if (outvar_on(14)) then

  do i = 1, gridcount
    rvar1d(i) = sv(i)%outvars%NPP
    if (.not.sv(i)%validcell) rvar1d(i) = rmissing
  end do

  call putvar1d(ofid,tpos,'NPP',gridstart,gridcount,rvar1d)

end if

!----------------------------
! Treecover
if (outvar_on(21)) then

  do i = 1, gridcount
    rvar1d(i) = sv(i)%outvars%treecover
    if (.not.sv(i)%validcell) rvar1d(i) = rmissing
  end do

  call putvar1d(ofid,tpos,'treecover',gridstart,gridcount,rvar1d)

end if

!----------------------------
! Grasscover
if (outvar_on(22)) then

  do i = 1, gridcount
    rvar1d(i) = sv(i)%outvars%grasscover
    if (.not.sv(i)%validcell) rvar1d(i) = rmissing
  end do

  call putvar1d(ofid,tpos,'grasscover',gridstart,gridcount,rvar1d)

end if


!----------------------------
! Cover
if (outvar_on(23)) then

  allocate(rvar2d(gridcount,npft))

  do i = 1, gridcount
    do j = 1, npft
      rvar2d(i,j) = sv(i)%outvars%cover(j)
      if (.not.sv(i)%validcell) rvar2d(i,:) = rmissing
    end do
  end do

  call putvar2d(ofid,tpos,'cover',gridstart,gridcount,rvar2d)

  deallocate(rvar2d)

end if

!----------------------------
! Nind
if (outvar_on(25)) then

  allocate(rvar2d(gridcount,npft))

  do i = 1, gridcount
    do j = 1, npft
      rvar2d(i,j) = sv(i)%outvars%nind(j)
      if (.not.sv(i)%validcell) rvar2d(i,:) = rmissing
    end do
  end do

  call putvar2d(ofid,tpos,'nind',gridstart,gridcount,rvar2d)

  deallocate(rvar2d)

end if

!----------------------------
! Crownarea
if (outvar_on(26)) then

  allocate(rvar2d(gridcount,npft))

  do i = 1, gridcount
    do j = 1, npft
      rvar2d(i,j) = sv(i)%outvars%crownarea(j)
      if (.not.sv(i)%validcell) rvar2d(i,:) = rmissing
    end do
  end do

  call putvar2d(ofid,tpos,'crownarea',gridstart,gridcount,rvar2d)

  deallocate(rvar2d)

end if

!----------------------------
! Soil moisture
if (outvar_on(27)) then

  allocate(rvar2d(gridcount,nl))

  do i = 1, gridcount
    do j = 1, nl
      rvar2d(i,j) = sv(i)%outvars%soilmoisture(j)
      if (.not.sv(i)%validcell) rvar2d(i,:) = rmissing
    end do
  end do

  call putvar2d(ofid,tpos,'soilmoisture',gridstart,gridcount,rvar2d)

  deallocate(rvar2d)

end if

!----------------------------
! Soil temperature
if (outvar_on(28)) then

  allocate(rvar2d(gridcount,nl))

  do i = 1, gridcount
    do j = 1, nl
      rvar2d(i,j) = sv(i)%outvars%soiltemp(j)
      if (.not.sv(i)%validcell) rvar2d(i,:) = rmissing
    end do
  end do

  call putvar2d(ofid,tpos,'soiltemp',gridstart,gridcount,rvar2d)

  deallocate(rvar2d)

end if

!----------------------------
! Height
if (outvar_on(46)) then

  allocate(rvar2d(gridcount,npft))

  do i = 1, gridcount
    do j = 1, npft
      rvar2d(i,j) = sv(i)%outvars%height(j)
      if (.not.sv(i)%validcell) rvar2d(i,:) = rmissing
    end do
  end do

  call putvar2d(ofid,tpos,'height',gridstart,gridcount,rvar2d)

  deallocate(rvar2d)

end if

!----------------------------
! GDD5
if (outvar_on(47)) then

  do i = 1, gridcount
    rvar1d(i) = sv(i)%outvars%GDD5
    if (.not.sv(i)%validcell) rvar1d(i) = rmissing
  end do

  call putvar1d(ofid,tpos,'GDD5',gridstart,gridcount,rvar1d)

end if

!--------------------------------------------
! flush the write buffer every 20 years of run

if (mod(tpos,20) == 0) ncstat = nf90_sync(ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)
!
! tpos = tpos + 1

deallocate(rvar1d)

end subroutine netcdf_output

!----------------------------------------------------------------------------------------------------------------------------

subroutine putvar1d(ofid,tpos,varname,gridstart,gridcount,rvar1d)

use parametersmod, only : i4,sp
use errormod,      only : ncstat,netcdf_err
use netcdf
use mpi

implicit none

integer(i4),                intent(in) :: ofid
integer(i4),                intent(in) :: tpos
character(*),               intent(in) :: varname
integer(i4),                intent(in) :: gridstart
integer(i4),                intent(in) :: gridcount
real(sp),     dimension(:), intent(in) :: rvar1d

! real(sp), allocatable, dimension(:,:) :: rvar1d

integer :: varid

!-----------------

!select type(var1d)
!type is (real)

! allocate(rvar1d(gridcount,1))

! rvar1d = reshape(var1d,[cntx,cnty])
!
! where (.not. cellmask) rvar2d = rmissing

! ncstat = nf90_open('./output/trash.nc',nf90_nowrite,ofid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,varname,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,rvar1d,start=[gridstart,tpos],count=[gridcount,1])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! deallocate(rvar2d)

!class default

!  write(stdout,*)'error an output data type was requested that is not supported'
!  stop

!end select

end subroutine putvar1d

!--------------------------------------------

subroutine putvar2d(ofid,tpos,varname,gridstart,gridcount,rvar2d)

use parametersmod, only : i4,sp
use errormod,      only : ncstat,netcdf_err
use netcdf

implicit none

integer(i4),                  intent(in) :: ofid
integer(i4),                  intent(in) :: tpos
character(*),                 intent(in) :: varname
integer(i4),                  intent(in) :: gridstart
integer(i4),                  intent(in) :: gridcount
real(sp),     dimension(:,:), intent(in) :: rvar2d

! real(sp), allocatable, dimension(:,:) :: rvar2d

integer :: l
integer :: varid

!-----------------

l = size(rvar2d,dim=2)

! allocate(rvar2d(gridcount,l))

ncstat = nf90_inq_varid(ofid,varname,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,rvar2d,start=[gridstart,1,tpos],count=[gridcount,l,1])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! deallocate(rvar2d)

end subroutine putvar2d

!--------------------------------------------

! subroutine putvar3d(ofid,tpos,varname,gridstart,gridcount,rvar1d)
!
! use parametersmod, only : i4,sp
! use errormod,      only : ncstat,netcdf_err
! use netcdf
!
! implicit none
!
! integer(i4),                intent(in) :: ofid
! integer(i4),                intent(in) :: tpos
! character(*),               intent(in) :: varname
! integer(i4),                intent(in) :: gridstart
! integer(i4),                intent(in) :: gridcount
! real(sp),     dimension(:), intent(in) :: rvar1d
!
! real(sp), allocatable, dimension(:,:) :: rvar2d
!
! integer :: varid
! integer :: i,l
! integer :: x,y
!
! !-----------------
!
! l = size(var2d,dim=2)
!
! !select type(var2d)
! !type is (real)
!
!   allocate(rvar3d(cntx,cnty,l))
!
!   y = 1
!   do i = 1,size(var2d,dim=1)
!
!     x = mod(i,cntx)
!     y = 1+ i / cntx
!
!     if (x==0) then
!       x = cntx
!       y = i/cntx
!     end if
!
!     rvar3d(x,y,:) = var2d(i,:)
!   end do
!
!   do i = 1,l
!     where (.not. cellmask) rvar3d(:,:,i) = rmissing
!   end do
!
!   ncstat = nf90_inq_varid(ofid,varname,varid)
!   if (ncstat/=nf90_noerr) call netcdf_err(ncstat)
!
!   ncstat = nf90_put_var(ofid,varid,rvar3d,start=[1,1,1,tpos],count=[cntx,cnty,l,1])
!   if (ncstat/=nf90_noerr) call netcdf_err(ncstat)
!
!   deallocate(rvar3d)
!
! !class default
!
! !  write(stdout,*)'error an output data type was requested that is not supported'
! !  stop
!
! !end select
!
! end subroutine putvar3d

!----------------------------------------------------------------------------------------------------------------------------



subroutine netcdf_open_mpi(file,fid)

use errormod, only : ncstat,netcdf_err
use netcdf
use mpi

implicit none

character(*), intent(in)  :: file
integer,      intent(out) :: fid

ncstat = nf90_open(file,nf90_write,fid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine netcdf_open_mpi

!----------------------------------------------------------------------------------------------------------------------------

subroutine netcdf_close(fid)

use errormod, only : ncstat,netcdf_err
use netcdf

implicit none

integer, intent(in) :: fid

ncstat = nf90_close(fid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

end subroutine netcdf_close

!----------------------------------------------------------------------------------------------------------------------------

end module netcdfoutputmod
