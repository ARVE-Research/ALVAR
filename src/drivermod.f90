module drivermod

! Module to recieve MPI job info, initiate states and distribute gridcell-level calculations

use parametersmod, only : i4,sp

implicit none

! Module variable for printgrid output
real(sp),      allocatable, dimension(:) :: plon
real(sp),      allocatable, dimension(:) :: plat
integer(i4),   allocatable, dimension(:) :: gruns
character(50), allocatable, dimension(:) :: textfilenames
integer(i4)                              :: ngrid
integer(i4)                              :: outmode
integer(i4)                              :: ofid
logical,                  dimension(100) :: outvar_on

contains

!-------------------------------------------------------

subroutine driver(info,job,rank)

use parametersmod,   only : i2,i4,sp,dp
use mpivarsmod,      only : mpivars
use pftparmod,       only : pftparameters
use printgridmod,    only : readgrdlist,setgrid,gentextfile
use randomdistmod,   only : randomstate
use statevarsmod,    only : inlon,inlat,time,indx,xlen,ylen,tlen,invars,nd,&
                            mon_metvars,day_metvars,soil_vars,gpp_vars,veg_vars,topo_vars,sv
use netcdfinputmod,  only : soildatainput,topodatainput
use netcdfoutputmod, only : netcdf_open_mpi,netcdf_close
use initmod,         only : initlonlat

implicit none

type(mpivars), target    , intent(in) :: info
integer(i4), dimension(2), intent(in) :: job
integer                  , intent(in) :: rank

! Pointers for mpi info variables
character(200), pointer :: outfile
character(200), pointer :: gridlistfile
character(200), pointer :: cfile_spinup
character(200), pointer :: cfile_transient
character(200), pointer :: soilfile
character(200), pointer :: topofile
logical,        pointer :: dospinup
logical,        pointer :: dotransient
integer(i4),    pointer :: spin_baseyr
integer(i4),    pointer :: spin_startyr
integer(i4),    pointer :: spinupyears
integer(i4),    pointer :: tran_baseyr
integer(i4),    pointer :: tran_startyr
integer(i4),    pointer :: transientyears
integer(i4),    pointer :: ilen
integer(i4),    pointer :: spin_tlen
integer(i4),    pointer :: tran_tlen

! Pointers for statevars
type(mon_metvars), pointer, dimension(:) :: monvars
type(day_metvars), pointer, dimension(:) :: dayvars
type(soil_vars),   pointer, dimension(:) :: soilvars
type(gpp_vars),    pointer, dimension(:) :: gppvars
type(veg_vars),    pointer, dimension(:) :: vegvars
type(topo_vars),   pointer, dimension(:) :: topovars
logical,           pointer, dimension(:) :: validcell
real(dp),          pointer, dimension(:) :: lon
real(dp),          pointer, dimension(:) :: lat
type(randomstate), pointer, dimension(:) :: georndst
real(sp),          pointer, dimension(:) :: elev

! Local variables
integer(i4), dimension(2) :: srt
integer(i4), dimension(2) :: cnt
character(200) :: infile
logical        :: spinup
integer(i4)    :: t0                ! start of time dimension for reading input file
integer(i4)    :: nt                ! count of time dimension for reading input file
integer(i4)    :: p0                ! time buffer-adjusted start position of input data array
integer(i4)    :: p1                ! time buffer-adjusted end position of input data array
integer(i4)    :: gridstart
integer(i4)    :: gridcount
integer(i4)    :: grid
integer(i4)    :: baseyr
integer(i4)    :: startyr
integer(i4)    :: calcyrs
integer(i4)    :: ndyear
integer(i4)    :: g
integer(i4)    :: yr
integer(i4)    :: d
integer(i4)    :: i
real(sp)       :: start_time
real(sp)       :: end_time

!-----------------------------------------------------------

call CPU_TIME(start_time)

!-------------------------

outfile         => info%outfile_list
gridlistfile    => info%gridlistfile
cfile_spinup    => info%cfile_spinup
cfile_transient => info%cfile_transient
soilfile        => info%soilfile
topofile        => info%topofile
dospinup        => info%dospinup
dotransient     => info%dotransient
spin_baseyr     => info%spin_baseyr
spin_startyr    => info%spin_startyr
spinupyears     => info%spinupyears
tran_baseyr     => info%tran_baseyr
tran_startyr    => info%tran_startyr
transientyears  => info%transientyears
ilen            => info%ilen
spin_tlen       => info%spin_tlen
tran_tlen       => info%tran_tlen

outmode         =  info%outmode
outvar_on       =  info%outvar_on

!-------------------------
! Create start and count array from the mpi job info

gridstart = job(1)
gridcount = job(2)

allocate(sv(gridcount))

monvars   => sv%monvars
soilvars  => sv%soilvars
topovars  => sv%topovars
vegvars   => sv%vegvars
lon       => sv%lon
lat       => sv%lat
georndst  => sv%georndst
elev      => sv%topovars%elev
validcell => sv%validcell

!-------------------------

call soildatainput(soilfile,gridstart,gridcount,soilvars)                             ! Read in soil variables to statevars

call topodatainput(topofile,gridstart,gridcount,topovars,validcell)                   ! Read in topographic variables to statevars

call pftparameters()                                                                  ! Initialize PFT-specific parameters

call initlonlat(cfile_spinup,gridstart,gridcount,xlen,ylen,inlon,inlat,indx,lon,lat)  ! Save lon and lat for each indexed grid

!-------------------------

if (outmode == 1) call readgrdlist(gridlistfile,plon,plat)

if (outmode == 1) call setgrid(lon,lat,plon,plat,ngrid,gruns)

! if (outmode == 1) call gentextfile(lon,lat,gruns,textfilenames)

if (outmode == 0) ngrid = gridcount

!-----------------------------------------------------------

infile  = cfile_spinup
baseyr  = spin_baseyr
startyr = spin_startyr
calcyrs = spinupyears
spinup  = .true.

if (dospinup) then

  call alvar(infile,outfile,rank,gridstart,gridcount,baseyr,startyr,calcyrs,spinup)

  deallocate(nd)
  deallocate(invars)

  spinup  = .false.

end if

!-----------------------------------------------------------

infile  = cfile_transient
baseyr  = tran_baseyr
startyr = tran_startyr
calcyrs = transientyears

if (dotransient) then

  call netcdf_open_mpi(outfile,ofid)

  call alvar(infile,outfile,rank,gridstart,gridcount,baseyr,startyr,calcyrs,spinup)

  call netcdf_close(ofid)

end if

!-----------------------------------------------------------

call CPU_TIME(end_time)

write(0,*) 'Rank:',rank, 'time spent on model:', end_time - start_time

!-------------------------

end subroutine driver

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------

subroutine alvar(infile,outfile,rank,gridstart,gridcount,baseyr,startyr,calcyrs,spinup)

use parametersmod,   only : i4,sp
use statevarsmod,    only : inlon,inlat,time,indx,xlen,ylen,invars,nd,sv
use orbitmod,        only : orbitpars,orbit,calcorbitpars
use netcdfinputmod,  only : metdatainput
use netcdfoutputmod, only : netcdf_output
use initmod,         only : initdate,initinvars,calcndyear
use printgridmod,    only : printgrid
use alvarmod,        only : alvar_annual
use alvarmod_daily,  only : alvar_daily

implicit none

character(200), intent(in) :: infile        ! Climate file input
character(200), intent(in) :: outfile       ! Output file
integer(i4),    intent(in) :: rank          ! Rank of mpi process
integer(i4),    intent(in) :: gridstart     ! Index of gridcell start
integer(i4),    intent(in) :: gridcount     ! Number of gridcell sent to process
integer(i4),    intent(in) :: baseyr        ! Base year of input climate file
integer(i4),    intent(in) :: startyr       ! Start year of run
integer(i4),    intent(in) :: calcyrs       ! Number of simulation years
logical,        intent(in) :: spinup        ! Spinup / transient logical

real(sp), pointer, dimension(:) :: elev

integer(i4), dimension(2) :: srt
integer(i4), dimension(2) :: cnt
integer(i4) :: tlen
integer(i4) :: t0
integer(i4) :: p0
integer(i4) :: p1
integer(i4) :: nt
integer(i4) :: cntt
integer(i4) :: grid
integer(i4) :: g
integer(i4) :: ndyear
integer(i4) :: yr
integer(i4) :: d
integer(i4) :: i

elev => sv%topovars%elev

call initdate(infile,rank,baseyr,startyr,calcyrs,tlen,t0,p0,p1,nt,cntt,nd)    ! Initialize the start, count and date of the job

srt = [gridstart, t0]    ! Start of gridcell, start of time dimension
cnt = [gridcount, nt]    ! Count of gridcell, count of time dimension (including +/- 1 buffer years)

!-------------------------

call initinvars(gridcount,cntt,invars)    ! Initilize the dimensions of the metvars variables (allocate by gridcount)

!-------------------------

call metdatainput(infile,gridcount,srt,cnt,cntt,p0,p1,xlen,ylen,tlen, &   ! Read in full array of monthly met input variable series
                  inlon,inlat,time,indx,elev,invars)

!-----------------------------------------------------------

yearloop : do yr = 1, calcyrs

  !-------------------------

  if (outmode == 0) print *, "Rank ", rank, " Working on year ", yr, ' out of ', calcyrs

  !-------------------------

  call calcndyear(yr,nd,ndyear)

  call calcorbitpars(startyr,yr,orbit)

  gridloop1 : do g = 1, ngrid

    if (outmode == 0) grid = g
    if (outmode == 1) grid = gruns(g)

    call alvar_annual(yr,grid,ndyear,calcyrs,spinup)

  end do gridloop1

  !-------------------------

  dayloop : do d = 1, ndyear

    diurnalloop : do i = 1, 2

      gridloop2 : do g = 1, ngrid

        if (outmode == 0) grid = g
        if (outmode == 1) grid = gruns(g)

        call alvar_daily(yr,grid,ndyear,calcyrs,d,i,spinup)

        ! if (outmode == 1) call printgrid(calcyrs,yr,d,i,grid,textfilenames(g),sv(grid))

      end do gridloop2

    end do diurnalloop

  end do dayloop

  !-------------------------
  ! Output variables into netcdf file in parallel
  if (outmode == 0 .and. .not.spinup) call netcdf_output(ofid,gridstart,gridcount,yr,outvar_on,sv)
  ! if (outmode == 0 .and. yr == calcyrs .and. .not.spinup) call netcdfoutput_onelayer(outfile,srt,cnt,yr,ndyear,nd,sv)

end do yearloop

end subroutine alvar

end module drivermod
