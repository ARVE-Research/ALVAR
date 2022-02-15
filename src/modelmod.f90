module modelmod

! Module to recieve MPI job info and distribute calculations to subroutines

use parametersmod,   only : i2,i4,sp,dp
use outputmod,       only : infompi,printgrid
use randomdistmod,   only : genrndstate
use pftparmod,       only : npft,pftparameters
use statevarsmod,    only : monvars,dayvars,soilvars,gppvars,vegvars,startyr,calcyrs,genvars,dayvars,topovars, &
                            ndyear,srt,cnt,clon,clat,gridlon,gridlat,lprint,gprint
use drivermod,       only : initdate,initlonlat,initmonvars,initsoilvars,initgppvars,initvegvars,inittopovars, &
                            initgeorndst,copygenvars,initdayvars,saveclonlat
use diurnaltempmod,  only : diurnaltemp,humidity
use orbitmod,        only : orbit,calcorbitpars
use radiationmod,    only : elev_Ratm,calcPjj,radpet,tdewpet,calcVPD,calctdew
use hourlyprecmod,   only : hourlyprec
use soilstatemod,    only : soilprep
use soilphysicsmod,  only : soilthermalprop,resistance,soiltemperature
use hydrologymod,    only : soilwater
use aetalphamod,     only : aet_alpha
use biome1mod,       only : initbiomevars,savebiomevars,calcbiome_year,calcbiome_mean
use phenologymod,    only : summerphenology
use waterbalancemod, only : simplesoilLPJ
use gppmod,          only : gpp
use nppmod,          only : npp
use establishmentmod,only : sapling,bioclim,establishment
use lightmod,        only : light
use turnovermod,     only : turnover
use killplantmod,    only : killplant
use allocationmod,   only : allocation
use mortalitymod,    only : mortality
use fireindexmod,    only : fireindex
use netcdfinputmod,  only : metdatainput,soildatainput,topodatainput,LAIdatainput
use netcdfoutputmod, only : netcdfoutput,netcdfoutput_onelayer
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

real(sp) :: start_time
real(sp) :: end_time

call CPU_TIME(start_time)

!-----------------------------------------------------------

call initdate(info,job,rank)        ! Initialize the start, count and date of the job

call initmonvars()                  ! Initilize the dimensions of the metvars variables (allocate by gridcount)

call initsoilvars()                 ! Initialize the dimensions of the soilvars variales (allocate by gridcount)

call initgppvars()                  ! Initialize the dimensions of the gppvars variables (allocate by gridcount and ndyear)

call initvegvars()                  ! Initialize the dimensions of the vegvars variales (allocate by gridcount)

call inittopovars()                 ! Initialize the dimensions of the topovar variales (allocate by gridcount)

call metdatainput(info)             ! Read in full array of monthly variable series

call soildatainput(info)            ! Read in full array of soil variables

call topodatainput(info)            ! Read in full array of topographic variables

call LAIdatainput(info)

call initlonlat(info,job,rank)      ! Save lon and lat for each indexed grid

call initgeorndst()                 ! Initialize dimension of the random state variables

call initbiomevars(calcyrs)         ! Initialize dimension of biomevars for saving annual met variables necessary for biome PFT (annual and long-term mean) estiamtes

call pftparameters()                ! Initialize PFT-specific parameters

!-----------------------------------------------------------

gridcount = job(2)

yearloop : do yr = 1, calcyrs

  ! Prepare soil physical properties for all grids from soildata input file
  gridloop_soilprep : do grid = 1, gridcount

    if (yr == 1) call soilprep(grid)

    if (yr == 1) call simplesoilLPJ(grid)       ! Currently stored in waterbalancemod.f90

  end do gridloop_soilprep

  !------

  ! Allocate dimension of dayvars from number of days in a year (365 or 366)
  call initdayvars(yr,gridcount)

  call calcorbitpars(startyr,yr,orbit)

  gridloop_metyear : do grid = 1, gridcount

    ! Generate an initial grid-specific (geohash) random state for weathergen subroutine
    ! 'georndst' variable in randomdistmod
    ! Only need to be called at year 1 as the grid random state will be saved after gwgen() as pass onto the next year
    if (yr == 1) call genrndstate(grid)

    ! Save current lon/lat of the grid in module variables 'clon' and 'clat' for later calculations
    ! NO LONGER IN USE (FOR NOW)
    ! call saveclonlat(grid)

    ! Copy 20 months (12 months +/- 4 months buffer) of monthly data for weathergen
    call copygenvars(yr,grid)

    ! Generate daily met variables from monthly series (of original input variables)
    call gwgen_new(grid)

  end do gridloop_metyear

  !------

  dayloop : do d = 1, ndyear

    gridloop_metday : do grid = 1, gridcount

      ! Copy 20 months (12 months +/- 4 months buffer) of monthly data
      call copygenvars(yr,grid)

      ! Subroutines to calculate other daily variables from tmin/tmax generated by gwgen()
      ! Calculate daylength, daytime and nighttime temperature
      call diurnaltemp(grid,d)

      ! Calculate atmospheric pressure variables based on elevation inputdata
      call elev_Ratm(grid,d)

      ! Calculate precipitation equitability index Pjj for tdew and pet in next subroutine
      call calcPjj(grid,d)

      ! Calculate shortwave (direct + diffuse) radiation and longwave radiation (stable solution loop with pet)
      ! Calculate potential evapotranspiration and dewpoint temperature
      call radpet(grid,d,1)

      ! call tdewpet(grid,d)         ! Routine written by Leo Lai (iterative stable solution tdew and PET routine, after Thronton et al., 2000; Kimball et al., 1997)

      ! Calculate relative humidity
      call humidity(grid,d)

      ! Calculate vapour pressure deficit
      call calcVPD(grid,d)

      ! Calculate actual evapotranspiration and AET/PET ration (alpha) using simple water bucket model (supply/demand from soil properties)
      call aet_alpha(yr,grid,d)

      ! Disaggregate 24 hour total prec into hourly series
      call hourlyprec(yr,grid,d)

      ! Calculate fire danger meter
      ! call fireindex(rank,yr,grid,d)

      ! Print grid data if the process recieved the user-specified lon/lat grid
      ! if (lprint .AND. grid == gprint) call printgrid(info,grid,yr,d)

    end do gridloop_metday

  end do dayloop

  !------

  ! Calculate biome type for first user-specified year, subroutine from BIOME1
  gridloop_biome : do grid = 1, gridcount

    call copygenvars(yr,grid)

    ! Save met variables for annual and long-term climatological mean biome PFT estimation
    call savebiomevars(yr,grid)

    ! Calculate the annual estimation of biome
    call calcbiome_year(yr,grid)

    if (yr == calcyrs) call calcbiome_mean(grid)      ! Calculate the long-term climate biome (only if arrived at last year loop)

  end do gridloop_biome

  !------

  gridloop_vegetation : do grid = 1, gridcount

    ! Copy 20 months (12 months +/- 4 months buffer) of monthly data
    call copygenvars(yr,grid)

    if (yr == 1) call sapling(yr,grid,1)

    call bioclim(yr,grid,1)

    call summerphenology(yr,grid,1)

    call establishment(yr,grid,1)

  end do gridloop_vegetation

  !------

  dayloop2 : do d = 1, ndyear

    diurnalloop : do i = 1, 2         ! 1 = day, 2 = night

      gridloop_soilwater : do grid = 1, gridcount

        ! if (soilvars(grid)%validcell) then

        ! call soilwater(rank,yr,grid,d,i)

        call soilthermalprop(grid)

        call resistance(grid,d)

        call soiltemperature(grid,d,i)

        if (i == 1) call gpp(yr,grid,d)

        if (i == 1) call npp(yr,grid,d)

      end do gridloop_soilwater

    end do diurnalloop

    if (d == ndyear) print *, rank,yr, d, &
                sum(gppvars%gpp_tot(1)), sum(gppvars%gpp_tot(2)), sum(gppvars%gpp_tot(3)),&
                sum(gppvars%gpp_tot(4)), sum(gppvars%gpp_tot(5)), &
                sum(gppvars%aresp(1)), sum(gppvars%aresp(2)), sum(gppvars%aresp(3)), &
                sum(gppvars%aresp(4)), sum(gppvars%aresp(5)), &
                sum(gppvars%npp_tot(1)), sum(gppvars%npp_tot(2)), sum(gppvars%npp_tot(3)), &
                sum(gppvars%npp_tot(4)), sum(gppvars%npp_tot(5))

  end do dayloop2


  !------

  ! dayloop3 : do d = 1, ndyear

  gridloop_vegetation2 : do grid = 1, gridcount

    ! Copy 20 months (12 months +/- 4 months buffer) of monthly data
    call copygenvars(yr,grid)

    call light(yr,grid,1)

    call turnover(yr,grid,1)

    ! call killplant(yr,grid,1)

    call allocation(yr,grid,ndyear)

    call light(yr,grid,1)

    ! call mortality(yr,grid,1)

  end do gridloop_vegetation2

  ! end do dayloop3

  !------

  ! Output variables into netcdf file in parallel
  ! call netcdfoutput(info,job,yr)
  if (yr == calcyrs) call netcdfoutput_onelayer(info,job,yr)

  deallocate(dayvars)

  if (rank == 0) write(0,*) 'Finished with year ', startyr+yr-1

end do yearloop

!-----------------------------------------------------------

call CPU_TIME(end_time)

write(0,*) 'Rank:',rank, 'time spent on model:', end_time - start_time


end subroutine model

end module modelmod
