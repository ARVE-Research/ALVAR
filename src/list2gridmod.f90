module list2gridmod

use parametersmod, only : i2,i4,sp

implicit none

integer(i2), parameter :: missing_i  = -32768
real(sp),    parameter :: missing_sp = -9999.

integer(i4), allocatable, dimension(:) :: xpos
integer(i4), allocatable, dimension(:) :: ypos

contains

!-------------------------------------------------------------------------------------

subroutine list2grid(info)

use parametersmod, only : i4,sp
use errormod,      only : ncstat,netcdf_err
use pftparmod,     only : npft
use mpivarsmod,    only : mpivars
use statevarsmod,  only : indx,xlen,ylen,nl
use netcdf

implicit none

type(mpivars), intent(in) :: info

logical, dimension(100) :: outvar_on

character(200) :: outfile           ! outfile name
character(200) :: outfile_list      ! outile list formmatted name

integer(i4) :: ifid
integer(i4) :: ofid
integer(i4) :: dimid
integer(i4) :: varid
! integer(i4) :: xlen
! integer(i4) :: ylen
integer(i4) :: ilen
integer(i4) :: tlen

integer, dimension(2) :: indx_loc

character(8) :: today
character(10) :: now

real(sp), allocatable, dimension(:) :: var_in_r
real(sp), allocatable, dimension(:,:) :: var_in_pft
real(sp), allocatable, dimension(:,:) :: var_out_r

real(sp), allocatable, dimension(:,:) :: var_out_ForFire

integer :: i
integer :: x, y

!-------------------------

ilen         = info%ilen
tlen         = info%transientyears
outfile      = info%outfile
outfile_list = info%outfile_list
outvar_on    = info%outvar_on

!-------------------------
! Save array for index xpos and ypos

allocate(xpos(ilen))
allocate(ypos(ilen))

do i = 1, ilen

  indx_loc = findloc(indx, i)

  xpos(i) = indx_loc(1)
  ypos(i) = indx_loc(2)

  ! print *, i, xpos(i), ypos(i)      ! Print to text file for saving

  if (mod(i,1000) == 0) write(0,*) 'Find loc at:', i, ' out of validcells ', ilen

end do

write(0,*)
write(0,*) 'xypos of all index found'

!-------------------------

ncstat = nf90_open(outfile_list,nf90_nowrite,ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_open(outfile,nf90_write,ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-------------------------

if (outvar_on(13)) call list2grid_2d(ifid,ofid,xlen,ylen,ilen,tlen,'livebiomass')

if (outvar_on(13)) call list2grid_2d(ifid,ofid,xlen,ylen,ilen,tlen,'GPP')

if (outvar_on(14)) call list2grid_2d(ifid,ofid,xlen,ylen,ilen,tlen,'NPP')

if (outvar_on(21)) call list2grid_2d(ifid,ofid,xlen,ylen,ilen,tlen,'treecover')

if (outvar_on(22)) call list2grid_2d(ifid,ofid,xlen,ylen,ilen,tlen,'grasscover')

if (outvar_on(23)) call list2grid_3d(ifid,ofid,xlen,ylen,ilen,tlen,npft,'cover')

if (outvar_on(25)) call list2grid_3d(ifid,ofid,xlen,ylen,ilen,tlen,npft,'nind')

if (outvar_on(26)) call list2grid_3d(ifid,ofid,xlen,ylen,ilen,tlen,npft,'crownarea')

if (outvar_on(27)) call list2grid_3d(ifid,ofid,xlen,ylen,ilen,tlen,nl,'soilmoisture')

if (outvar_on(28)) call list2grid_3d(ifid,ofid,xlen,ylen,ilen,tlen,nl,'soiltemp')

if (outvar_on(46)) call list2grid_3d(ifid,ofid,xlen,ylen,ilen,tlen,npft,'height')

if (outvar_on(47)) call list2grid_2d(ifid,ofid,xlen,ylen,ilen,tlen,'GDD5')

!-------------------------

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

write(0,*)
write(0,*) 'Finished all'
write(0,*)

!-------------------------

end subroutine list2grid

!-------------------------------------------------------------------------------------

subroutine list2grid_2d(ifid,ofid,xlen,ylen,ilen,tlen,varname)

use errormod,      only : ncstat,netcdf_err
use netcdf

implicit none

integer(i4),  intent(in) :: ifid
integer(i4),  intent(in) :: ofid
integer(i4),  intent(in) :: xlen
integer(i4),  intent(in) :: ylen
integer(i4),  intent(in) :: ilen
integer(i4),  intent(in) :: tlen
character(*), intent(in) :: varname


real(sp), allocatable, dimension(:,:)   :: var_in
real(sp), allocatable, dimension(:,:,:) :: var_out
integer(i4) :: varid
integer(i4) :: x
integer(i4) :: y
integer(i4) :: i

allocate(var_in(ilen,tlen))
allocate(var_out(xlen,ylen,tlen))

var_out = missing_sp

ncstat = nf90_inq_varid(ifid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do i = 1, ilen

  x = xpos(i)
  y = ypos(i)

  var_out(x,y,:) = var_in(i,:)

end do

!-------------------------

write(0,*)
write(0,*) 'Putting ', varname, ' into new output file'
write(0,*)

ncstat = nf90_put_var(ofid,varid,var_out)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!-------------------------

deallocate(var_in)
deallocate(var_out)

!-------------------------

end subroutine list2grid_2d

!-------------------------------------------------------------------------------------

subroutine list2grid_3d(ifid,ofid,xlen,ylen,ilen,tlen,dim3,varname)

use errormod,      only : ncstat,netcdf_err
use netcdf

implicit none

integer(i4),  intent(in) :: ifid
integer(i4),  intent(in) :: ofid
integer(i4),  intent(in) :: xlen
integer(i4),  intent(in) :: ylen
integer(i4),  intent(in) :: ilen
integer(i4),  intent(in) :: tlen
integer(i4),  intent(in) :: dim3
character(*), intent(in) :: varname


real(sp), allocatable, dimension(:,:,:)   :: var_in
real(sp), allocatable, dimension(:,:,:,:) :: var_out
integer(i4) :: varid
integer(i4) :: x
integer(i4) :: y
integer(i4) :: i

allocate(var_in(ilen,dim3,tlen))
allocate(var_out(xlen,ylen,dim3,tlen))

var_out = missing_sp

ncstat = nf90_inq_varid(ifid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do i = 1, ilen

  x = xpos(i)
  y = ypos(i)

  var_out(x,y,:,:) = var_in(i,:,:)

end do

!-------------------------

write(0,*)
write(0,*) 'Putting ', varname, ' into new output file'
write(0,*)

ncstat = nf90_put_var(ofid,varid,var_out)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!-------------------------

deallocate(var_in)
deallocate(var_out)

!-------------------------

end subroutine list2grid_3d

!-------------------------------------------------------------------------------------

end module list2gridmod
