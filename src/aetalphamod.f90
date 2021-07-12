module aetalphamod

use parametersmod, only : i4,sp,dp
use metvarsmod,    only : dayvars,cnt
use simplesoilmod, only : soilvars

implicit none

! Global variables
! integer(i4) :: gridcount = cnt(1)

real(sp), allocatable, dimension(:) :: m       		! Instantaneous soil moisture (mm)

public :: aet_alpha

contains

!-----------------------------------------

subroutine aet_alpha(year,grid,day)

implicit none

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day


!arguments

real(sp) :: awc             ! Total water holding capacity in the top meter of soil (mm); actually, total whc (mm3/mm3) multiplied by soil depth (mm) PMC 25 jan '10
real(sp), pointer :: dayl   ! Day length (hr)
real(sp), pointer :: prec   ! Daily precipitation (mm)
real(sp), pointer :: dpet   ! Daily potential evapotranspiration (mm)
real(sp), pointer :: daet   ! Daily actual evapotranspiration (mm)
real(sp), pointer :: alpha  ! Ratio of AET/PET (fraction)

!parameters

real(sp), parameter :: etmax = 1.0_sp  ! Maximum equilibrium evapotranspiration rate (mm h-1)

!local variables
real(sp) :: Emd     ! Emax: maximum evapotranspiration rate from saturated soils(mm d-1)
real(sp) :: supply  ! Water supply (mm)
real(sp) :: demand  ! Water demand (mm)

integer(i4) :: gridcount

!-------------------

awc   =  sum(soilvars(grid)%whc(1:5))    ! Top meter (i.e. first five layer) of water carrying capcity
dayl  => dayvars(grid,day)%dayl
prec  => dayvars(grid,day)%prec
dpet  => dayvars(grid,day)%dpet
daet  => dayvars(grid,day)%daet
alpha => dayvars(grid,day)%alpha

!------

! Initialize instatenous soil moisture for each gridcell
if (year == 1 .and. grid == 1 .and. day == 1) then

  gridcount = cnt(1)

  allocate(m(gridcount))

end if

if (year == 1 .and. day == 1) m(grid) = 0.

!------

! If loop to filter out missing soil data
if (soilvars(grid)%whc(1) == -9999.) then

  daet = -9999.
  alpha = -9999.

else

  ! Demand = PET (simple case assumption)
  demand = dpet

  ! First, calculate supply based on previous day's soil moisture
  Emd = etmax * dayl !first find total etmax for the whole day

  supply = min(Emd, m(grid)) !Emd * m / awc

  ! Second, calculate aet for this day using supply
  daet = min(supply, demand)

  if (dpet > 0.0) then

    alpha = daet / dpet

  else

    alpha = 1.0

  end if

  ! Update soil moisture to reflect current day's aet --> bucket soil moisture (Leo)
  m(grid) = min(m(grid)+(prec-daet), awc)

  !write(*,'(7f10.3)')prec,emd,m,supply,pet,aet,alpha

end if


end subroutine aet_alpha

!-----------------------------------------

end module aetalphamod
