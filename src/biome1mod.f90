module biome1mod

use parametersmod, only : i2,i4,sp,dp,missing_i2
use statevarsmod,    only : genvars,dayvars,ndyear,cnt,calcyrs

implicit none

public :: initbiomevars
public :: savebiomevars
public :: calcbiome_year
public :: calcbiome_mean

integer(i4), parameter :: npfts = 13

!---------------------------------------------------------------------

type biomedata

  real(dp) :: tcm       ! Temperature of coldest month in year (degC)
  real(dp) :: twm       ! Temperature of warmest month in year (degC)
  real(dp) :: gdd5      ! Growing degree days above 5 degC
  real(dp) :: gdd0      ! Growing degree days above 0 degC
  real(dp) :: alpha     ! Annual mean of alpha (AET/PET ratio)
  real(dp) :: aprec     ! Annual sum of precipitation (mm)
  real(dp) :: apet      ! Annual sum potential evapotranspiration (mm)

  integer(i4) :: ndyear   ! Number of days in the year

  integer(i4) :: biome    ! Biome PFT types (1-17)

end type biomedata

! Global variable to store biome variables for mean climate estimate of biome PFTs
! Allocate variable with gridcount and calcyears for calculations of mean
type(biomedata), allocatable, dimension(:,:) :: biomevars

integer(i4), allocatable, dimension(:) :: biome_mean

!---------------------------------------------------------------------

type pftparameters

  character(80) :: name
  real(dp) :: tcmax
  real(dp) :: tcmin
  real(dp) :: gdd5min
  real(dp) :: gdd0min
  real(dp) :: twmin
  real(dp) :: almin
  real(dp) :: almax

end type pftparameters

type(pftparameters), dimension(npfts) :: pftpar


contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine initbiomevars(calcyrs)

! Subroutine to allocate biomevars

implicit none

integer(i4), intent(in) :: calcyrs

integer(i4) :: gridcount

!------

gridcount = cnt(1)

allocate(biomevars(gridcount,calcyrs))

!------

end subroutine initbiomevars

!---------------------------------------------------------------------

subroutine setpftpars()

! Subroutine to set the values for pftpar global variable
! To be called in calcbiome subroutine if gird number = 1 (i.e., first grid in rank)

implicit none

pftpar(1)%name  = "warm broadleaf evergreen"
pftpar(2)%name  = "warm broadleaf deciduous"
pftpar(3)%name  = "temperate broadleaf evergreen"
pftpar(4)%name  = "temperate broadleaf deciduous"
pftpar(5)%name  = "temperate needleleaf evergreen"
pftpar(6)%name  = "cold needleleaf evergreen"
pftpar(7)%name  = "cold needle- and broadleaf deciduous"
pftpar(8)%name  = "sclerophyll/succulent"
pftpar(9)%name  = "warm herbaceous"
pftpar(10)%name = "cool herbaceous"
pftpar(11)%name = "cold herbaceous"
pftpar(12)%name = "hot shrub"
pftpar(13)%name = "cold shrub"

pftpar%tcmin   = [15.5,     15.5,     6.,   -15.,   -19.,   -35., -9999.,     5., -9999., -9999., -9999., -9999., -9999.]

pftpar%tcmax   = [-9999., -9999., -9999.,   15.5,     5.,    -2.,     5., -9999., -9999., -9999., -9999., -9999., -9999.]

pftpar%gdd5min = [-9999., -9999., -9999.,  1200.,   900.,   350.,   350., -9999., -9999.,   500., -9999., -9999., -9999.]

pftpar%gdd0min = [-9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999.,   100., -9999.,   100.]

pftpar%twmin   = [-9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999.,    22., -9999., -9999.,    22., -9999.]

pftpar%almin   = [0.80,     0.45,   0.65,   0.65,   0.65,   0.75,   0.65,   0.28,   0.18,   0.33,   0.33, -9999., -9999.]

pftpar%almax   = [-9999.,   0.95, -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999.]


end subroutine setpftpars

!---------------------------------------------------------------------

subroutine calcbiome_year(year,grid)

! Subroutine to calculate biome type based on met variables (code copied from BIOME1) - Leo O Lai (June 2021)
! Input variables modified for the ALVAR model

use statevarsmod, only : vegvars

implicit none

! Arguments
integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid

integer(i4), dimension(2) :: biome

real(dp) :: tcm       ! Temperature of coldest month in year (degC)
real(dp) :: twm       ! Temperature of warmest month in year (degC)
real(dp) :: gdd5      ! Growing degree days above 5 degC
real(dp) :: gdd0      ! Growing degree days above 0 degC
real(dp) :: alpha     ! Annual mean of alpha (AET/PET ratio)
real(dp) :: aprec     ! Annual sum of precipitation (mm)
real(dp) :: apet      ! Annual sum potential evapotranspiration (mm)

! Parameters
real(dp),    parameter :: ud = -9999.
integer(i4), parameter :: tests = 7

! Local variables
integer(i4) :: pft

integer(i4) :: gridcount
integer(i4) :: bit
integer(i4) :: bval

logical, dimension(npfts) :: present
logical, dimension(tests) :: passed

!------

if (year == 1 .and. grid == 1) call setpftpars()

!------

tcm   = biomevars(grid,year)%tcm
twm   = biomevars(grid,year)%twm
gdd5  = biomevars(grid,year)%gdd5
gdd0  = biomevars(grid,year)%gdd0
alpha = biomevars(grid,year)%alpha
aprec = biomevars(grid,year)%aprec
apet  = biomevars(grid,year)%apet

!------

biome = -9999.

!calculate present pfts

present = .true.

do pft = 1,npfts

  passed = .true.

  !tcmin
  if (pftpar(pft)%tcmin /= ud .and. tcm < pftpar(pft)%tcmin)      passed(1) = .false.

  !tcmax
  if (pftpar(pft)%tcmax /= ud .and. tcm > pftpar(pft)%tcmax)      passed(2) = .false.

  !GDD5
  if (pftpar(pft)%gdd5min /= ud .and. gdd5 < pftpar(pft)%gdd5min) passed(3) = .false.

  !GDD0
  if (pftpar(pft)%gdd0min /= ud .and. gdd0 < pftpar(pft)%gdd0min) passed(4) = .false.

  !twmin
  if (pftpar(pft)%twmin /= ud .and. twm < pftpar(pft)%twmin)      passed(5) = .false.

  !alpha min
  if (pftpar(pft)%almin /= ud .and. alpha < pftpar(pft)%almin)    passed(6) = .false.

  !alpha max
  if (pftpar(pft)%almax /= ud .and. alpha > pftpar(pft)%almax)    passed(7) = .false.

  if (all(passed)) then
    present(pft) = .true.
  else
    present(pft) = .false.
  end if

end do

!------

if (aprec - apet < -100.) then
  present(1:7) = .false.
end if

!------
! Choose biome based on present pfts
! Calculate an unique integer value that represents the PFTs present

bval = 0

do pft = 1,npfts

  if (present(pft)) then
    bit = 1
  else
    bit = 0
  end if

  bval = bval + bit * 2**(npfts-pft)

end do

!------

do

  if (present(1)) then
    if (present(2)) then
      biome = 2  !tropical raingreen
      exit
    else
      biome = 1  !tropical evergreen
      exit
    end if
  end if

  if (present(2)) then
    biome = 3   !tropical deciduous
    exit
  end if

  if (present(3)) then
    if (present(4)) then
      biome = 5   !temperate deciduous because temperate deciduous will always dominate in situations where it is present
      exit
    else
      biome = 4   !temperate evergreen (warm mixed)
      exit
    end if
  end if

  if (present(4)) then
    if (present(5) .and. present(6) .and. present(7)) then
      biome = 6  !cool mixed
      exit
    !else if (present(5) .and. present(7)) then
    else if (present(5) .and. tcm < 1.) then
      biome = 6
      exit
    else
      biome = 5  !temperate deciduous
      exit
    end if
  end if

  if (present(5) .and. present(6) .and. present (7)) then
    biome = 7   !cool conifer
    exit
  end if

  if (present(6) .and. present(7)) then
    biome = 8  !cold evergreen
    exit
  end if

  if (present(5) .and. present(7)) then
    biome = 9 !cold mixed
    exit
  end if

  if (present(7)) then
    biome = 10  !cold deciduous
    exit
  end if

  if (present(8)) then
    biome = 11  !xerophytic (sclerophyll)
    exit
  end if

  if (present(9)) then
    biome = 12  !warm grass
    exit
  end if

  if (present(10) .and. present(11)) then
    biome = 13  !cool grass
    exit
  else if (present(11)) then
    biome = 14  !tundra
    exit
  end if

  if (present(12)) then
    biome = 15  !hot desert
    exit
  end if

  if (present(13)) then
    biome = 16  !cool desert
    exit
  end if

  biome = 17  !barren
  exit

end do

!---

! write(0,*)tcm,gdd5,gdd0,twm,alpha
!write(0,*)present
! write(0,*) biome

! biome = bval

if (dayvars(grid,1)%daet == -9999.) biome = missing_i2

biomevars(grid,year)%biome = biome(1)
vegvars(grid)%biome = biome(1)


end subroutine calcbiome_year

!---------------------------------------------------------------------

subroutine calcbiome_mean(grid)

! Subroutine to calculate biome type based on met variables (code copied from BIOME1) - Leo O Lai (June 2021)
! Input variables modified for the ALVAR model

implicit none

integer(i4), intent(in) :: grid

integer(i4), dimension(2) :: biome

real(dp) :: tcm       ! Temperature of coldest month in year (degC)
real(dp) :: twm       ! Temperature of warmest month in year (degC)
real(dp) :: gdd5      ! Growing degree days above 5 degC
real(dp) :: gdd0      ! Growing degree days above 0 degC
real(dp) :: alpha     ! Annual mean of alpha (AET/PET ratio)
real(dp) :: aprec     ! Annual sum of precipitation (mm)
real(dp) :: apet      ! Annual sum potential evapotranspiration (mm)

! Parameters
real(dp),    parameter :: ud = -9999.
integer(i4), parameter :: tests = 7

! Local variables
integer(i4) :: pft

integer(i4) :: gridcount
integer(i4) :: bit
integer(i4) :: bval

logical, dimension(npfts) :: present
logical, dimension(tests) :: passed

integer :: i

!------

if (grid == 1) then

  call setpftpars()

  gridcount = cnt(1)

  allocate(biome_mean(gridcount))

end if

!------

tcm   = sum(biomevars(grid,:)%tcm) / calcyrs
twm   = sum(biomevars(grid,:)%twm)  / calcyrs
gdd5  = sum(biomevars(grid,:)%gdd5) / calcyrs
gdd0  = sum(biomevars(grid,:)%gdd0) / calcyrs
aprec = sum(biomevars(grid,:)%aprec) / calcyrs
apet  = sum(biomevars(grid,:)%apet) / calcyrs

alpha = 0.

do i = 1, calcyrs

  alpha = alpha + biomevars(grid,i)%alpha * biomevars(grid,i)%ndyear

end do

alpha = alpha / sum(biomevars(grid,:)%ndyear)

!------

biome = -9999.

!calculate present pfts

present = .true.

do pft = 1,npfts

  passed = .true.

  !tcmin
  if (pftpar(pft)%tcmin /= ud .and. tcm < pftpar(pft)%tcmin)      passed(1) = .false.

  !tcmax
  if (pftpar(pft)%tcmax /= ud .and. tcm > pftpar(pft)%tcmax)      passed(2) = .false.

  !GDD5
  if (pftpar(pft)%gdd5min /= ud .and. gdd5 < pftpar(pft)%gdd5min) passed(3) = .false.

  !GDD0
  if (pftpar(pft)%gdd0min /= ud .and. gdd0 < pftpar(pft)%gdd0min) passed(4) = .false.

  !twmin
  if (pftpar(pft)%twmin /= ud .and. twm < pftpar(pft)%twmin)      passed(5) = .false.

  !alpha min
  if (pftpar(pft)%almin /= ud .and. alpha < pftpar(pft)%almin)    passed(6) = .false.

  !alpha max
  if (pftpar(pft)%almax /= ud .and. alpha > pftpar(pft)%almax)    passed(7) = .false.

  if (all(passed)) then
    present(pft) = .true.
  else
    present(pft) = .false.
  end if

end do

!------

if (aprec - apet < -100.) then
  present(1:7) = .false.
end if

!------
! Choose biome based on present pfts
! Calculate an unique integer value that represents the PFTs present

bval = 0

do pft = 1,npfts

  if (present(pft)) then
    bit = 1
  else
    bit = 0
  end if

  bval = bval + bit * 2**(npfts-pft)

end do

!------

do

  if (present(1)) then
    if (present(2)) then
      biome = 2  !tropical raingreen
      exit
    else
      biome = 1  !tropical evergreen
      exit
    end if
  end if

  if (present(2)) then
    biome = 3   !tropical deciduous
    exit
  end if

  if (present(3)) then
    if (present(4)) then
      biome = 5   !temperate deciduous because temperate deciduous will always dominate in situations where it is present
      exit
    else
      biome = 4   !temperate evergreen (warm mixed)
      exit
    end if
  end if

  if (present(4)) then
    if (present(5) .and. present(6) .and. present(7)) then
      biome = 6  !cool mixed
      exit
    !else if (present(5) .and. present(7)) then
    else if (present(5) .and. tcm < 1.) then
      biome = 6
      exit
    else
      biome = 5  !temperate deciduous
      exit
    end if
  end if

  if (present(5) .and. present(6) .and. present (7)) then
    biome = 7   !cool conifer
    exit
  end if

  if (present(6) .and. present(7)) then
    biome = 8  !cold evergreen
    exit
  end if

  if (present(5) .and. present(7)) then
    biome = 9 !cold mixed
    exit
  end if

  if (present(7)) then
    biome = 10  !cold deciduous
    exit
  end if

  if (present(8)) then
    biome = 11  !xerophytic (sclerophyll)
    exit
  end if

  if (present(9)) then
    biome = 12  !warm grass
    exit
  end if

  if (present(10) .and. present(11)) then
    biome = 13  !cool grass
    exit
  else if (present(11)) then
    biome = 14  !tundra
    exit
  end if

  if (present(12)) then
    biome = 15  !hot desert
    exit
  end if

  if (present(13)) then
    biome = 16  !cool desert
    exit
  end if

  biome = 17  !barren
  exit

end do

!---

! write(0,*)tcm,gdd5,gdd0,twm,alpha
!write(0,*)present
! write(0,*) biome

! biome = bval

if (dayvars(grid,1)%daet == -9999.) biome = missing_i2

biome_mean(grid) = biome(1)


end subroutine calcbiome_mean

!---------------------------------------------------------------------

subroutine savebiomevars(year,grid)

! Subroutine to save annual met variables into biomevars for mean climate estimate of
! biome PFTs

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid

biomevars(grid,year)%tcm = minval(genvars%tmp(5:16))
biomevars(grid,year)%twm = maxval(genvars%tmp(5:16))

biomevars(grid,year)%gdd5  = sum(dayvars(grid,1:ndyear)%tmean - 5.0, mask = dayvars(grid,1:ndyear)%tmean > 5.0)
biomevars(grid,year)%gdd0  = sum(dayvars(grid,1:ndyear)%tmean - 0.0, mask = dayvars(grid,1:ndyear)%tmean > 0.0)

biomevars(grid,year)%alpha = sum(dayvars(grid,1:ndyear)%alpha) / ndyear

biomevars(grid,year)%aprec = sum(dayvars(grid,1:ndyear)%prec)
biomevars(grid,year)%apet  = sum(dayvars(grid,1:ndyear)%dpet)

biomevars(grid,year)%ndyear = ndyear


end subroutine savebiomevars

!---------------------------------------------------------------------

end module biome1mod
