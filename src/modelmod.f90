module modelmod

! Module to recieve MPI job info and distribute calculations to subroutines

use parametersmod,   only : i2,i4,sp,dp
use outputmod,       only : infompi,printgrid
use randomdistmod,   only : genrndstate
use metvarsmod,      only : monvars,dayvars,startyr,calcyrs,genvars,dayvars,ndyear,srt,cnt,clon,clat,gridlon,gridlat,lprint,gprint
use drivermod,       only : initdate,initlonlat,initmonvars,copygenvars,initdayvars,saveclonlat
use diurnaltempmod,  only : diurnaltemp,humidity,calctdew
use orbitmod,        only : orbit,calcorbitpars
use radiationmod,    only : calcPjj,radpet,tdewpet,calcVPD
use netcdfinputmod,  only : netcdfinput
use netcdfoutputmod, only : netcdfoutput
use gwgenmod,        only : gwgen
use gwgenmodnew,     only : gwgen_new
use netcdf
use mpi

implicit none

contains

!-------------------------------------------------------
subroutine model(info,job,rank)

type(infompi), target    , intent(in) :: info
integer(i4), dimension(2), intent(in) :: job
integer                  , intent(in) :: rank

! Pointers for mpi info variables
character(100), pointer :: infile
character(100), pointer :: outfile
character(100), pointer :: timestring
integer(i4)   , pointer :: nproc
integer(i4)   , pointer :: t0
integer(i4)   , pointer :: nt

integer(i4) :: gridcount
integer :: grid
integer :: yr
integer :: d
integer :: i

real :: start_time, end_time
integer :: test

call CPU_TIME(start_time)

!-----------------------------------------------------------

call initdate(info,job,rank)        ! Initialize the start, count and date of the job

call initmonvars()                  ! Initilize the dimensions of the metvars variables

call netcdfinput(info)              ! Read in full array of monthly variable series

call initlonlat(info,job,rank)

!-----------------------------------------------------------

gridcount = job(2)

yearloop : do yr = 1, calcyrs

  ! Allocate dimension of dayvars from number of days in a year (365 or 366)
  call initdayvars(yr,gridcount)

  gridloop : do grid = 1, gridcount

    ! Generate an initial grid-specific (geohash) random state for weathergen subroutine
    ! 'georndst' variable in randomdistmod
    call genrndstate(grid)

    ! Save current lon/lat of the grid in module variables 'clon' and 'clat' for later calculations
    ! call saveclonlat(grid)

    ! Copy 20 months (12 months +/- 4 months buffer) of monthly data for weathergen
    call copygenvars(yr,grid)

    ! Generate daily met variables from monthly series (of original input variables)
    call gwgen_new(grid)

    ! do d = 1, ndyear
    !   if (clon == -102.25 .AND. clat == 37.75) then
    !     print *,yr,d, clon,clat, dayvars(grid,d)%tmin,dayvars(grid,d)%tmax, dayvars(grid,d)%cldf
    !   end if
    ! end do

  end do gridloop

  !--------

  dayloop : do d = 1, ndyear

    gridloop2 : do grid = 1, gridcount

      call diurnaltemp(grid,d)

      call calcPjj(grid,d)

      call radpet(grid,d,1.,1)

      call calctdew(grid,d)        ! Routine written by Leo Lai (after Kimbell et al., 1997)

      ! call tdewpet(grid,d)       ! Routine written by Leo Lai (iterative tdew and PET routine, after Thronton et al., 2000)

      call humidity(grid,d)

      call calcVPD(grid,d)

      ! Print grid data if the process recieved the user-specified lon/lat grid
      if (lprint .AND. grid == gprint) call printgrid(info,grid,yr,d)

    end do gridloop2

  end do dayloop

  ! Output dayvars into netcdf file in parallel
  call netcdfoutput(info,job,yr)

  deallocate(dayvars)

end do yearloop







!-----------------------------------------------------------

call CPU_TIME(end_time)

write(0,*) 'Rank:',rank, 'time spent on model:', end_time - start_time

! write(0,*) 'Test:', test


end subroutine model

end module modelmod
