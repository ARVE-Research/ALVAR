module pftparmod

use parametersmod, only : i4,sp,dp

implicit none

integer, parameter :: npft = 5
integer, parameter :: tnpft = 5       ! Temporary quick fix variable for this module in case for testing out PFT combinations

real(sp), dimension(30,tnpft) :: pftpar

logical, dimension(tnpft) :: tree
logical, dimension(tnpft) :: evergreen
logical, dimension(tnpft) :: summergreen
logical, dimension(tnpft) :: raingreen
logical, dimension(tnpft) :: needle
logical, dimension(tnpft) :: boreal
logical, dimension(tnpft) :: c4

contains

subroutine pftparameters()

! Subroutine to initilize all PFT parameters

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

! 1) fraction of roots in upper soil layer --- rootprop called in gppmod
pftpar(1,:) = [0.85, 0.80, 0.90, 0.90, 0.90]

! 2) C3 or C4 photosynthesis pathway, C4 = 1.0
pftpar(2,:) = [0.0, 0.0, 0.0, 0.0, 1.0]

! 3) water scale value at which leaves shed by drought deciduous PFT --- called in gppmod
pftpar(3,:) = [0.0, 0.0, 0.0, 0.35, 0.35]

! 4) canopy conductance component (gmin, mm/s) not associated with photosynthesis (Haxeltine & Prentice 1996, Table 4) --- called in gppmod
pftpar(4,:) = [0.50, 0.50, 0.30, 0.50, 0.50]

! 5) maintenance respiration coefficient --- called in nppmod
pftpar(5,:) = [0.10, 1.00, 1.20, 0.70, 0.20]

! 6) maximum foliar N content (mg/g) --- called in gppmod
pftpar(6,:) = [100.0, 120.0, 100.0, 100.0, 100.0]

! 7) leaf longevity (yr) --- for calculating SLA
pftpar(7,:) = [2.0, 0.5, 2.0, 1.0, 1.0]

! 8) leaf turnover period (yr) --- called in turnovermod
pftpar(8,:) = [2.0, 1.0, 2.0, 1.0, 1.0]

! 9) sapwood turnover period to heartwood (yr) --- called in turnovermod
pftpar(9,:) = [20.0, 20.0, 20.0, 1.0, 1.0]

! 10) root turnover period (yr) --- called in turnovermod
pftpar(10,:) = [2.0, 1.0, 2.0, 2.0, 2.0]

! 11) leaf C:N mass ratio --- called in nppmod
pftpar(11,:) = [29.0, 29.0, 29.0, 29.0, 29.0]

! 12) sapwood C:N mass ratio --- called in nppmod
pftpar(12,:) = [330.0, 330.0, 330.0, 0.0, 0.0]

! 13) root C:N mass ratio --- called in nppmod
pftpar(13,:) = [29.0, 29.0, 29.0, 29.0, 29.0]

! 14) leaf type: broadleaved (1), needleleaved (2) or grass (3)
pftpar(14,:) = [1.0, 1.0, 2.0, 3.0, 3.0]

! 15) phenology type: evergreen (1), summergreen (2), raingreen (3), any type (4)
pftpar(15,:) = [1.0, 2.0, 1.0, 4.0, 4.0]

! 16) leaf to root ratio under non-water stressed conditions --- called in allocationmod
pftpar(16,:) = [1.0, 1.0, 1.0, 0.75, 0.75]

! 17) summergreen phenology ramp, GDD5 requirement to grow full leaf canopy --- called in phenologymod
pftpar(17,:) = [1000.0, 200.0, 1000.0, 100.0, 100.0]

! 18) tree maximum crown area (m2)
pftpar(18,:) = [30.0, 30.0, 30.0, 0.0, 0.0]

! 19) sapling / grass initialization LAI
pftpar(19,:) = [4.0, 4.0, 4.0, 0.001, 0.001]

! 20) sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)] ratio
pftpar(20,:) = [3.0, 3.0, 3.0, 3.0, 3.0]

! 21) boreal PFT (1), non-boreal PFT (0)
pftpar(21,:) = [0.0, 0.0, 1.0, 1.0, 0.0]

! 22) low temperature limit for CO2 update --- called in gppmod
pftpar(22,:) = [2.0, -4.0, -4.0, -4.0, 6.0]

! 23) lower range of temperature optimum for photosynthesise --- called in gppmod
pftpar(23,:) = [25.0, 20.0, 15.0, 10.0, 20.0]

! 24) upper range of temperature optimum for photosynthesis --- called in gppmod
pftpar(24,:) = [30.0, 25.0, 25.0, 30.0, 45.0]

! 25) high temperature limit for CO2 update --- called in gppmod
pftpar(25,:) = [55.0, 38.0, 38.0, 45.0, 55.0]

! 26) optimal Ci/Ca ratio --- lambda called in gppmod
pftpar(26,:) = [0.90, 0.80, 0.80, 0.65, 0.40]

! 27) minimum coldest monthly mean temperature --- called in bioclim
pftpar(27,:) = [15.5, -17.0, -32.5, -1000.0, 15.5]

! 28) maximum coldest monthly mean temperature --- called in bioclim
pftpar(28,:) = [1000.0, 15.0, -2.0, 15.5, 1000.0]

! 29) minimum growing degree days (at or above 5 deg C) --- called in bioclim
pftpar(29,:) = [0.0, 1200.0, 600.0, 0.0, 0.0]

! 30) upper limit of temperature of the warmest month
pftpar(30,:) = [1000.0, 1000.0, 23.0, 1000.0, 1000.0]

!---------------------------------------------------------------------

tree        = [ .true.,  .true.,  .true.,  .false., .false. ]
evergreen   = [ .true.,  .false., .true.,  .false., .false. ]
summergreen = [ .false., .true.,  .false., .false., .false. ]
raingreen   = [ .false., .false., .false., .false., .false. ]
needle      = [ .false., .false., .true.,  .false., .false. ]
boreal      = [ .false., .false., .true.,  .true.,  .false. ]
c4          = [ .false., .false., .false., .false., .true.  ]

!---------------------------------------------------------------------

end subroutine pftparameters



end module pftparmod



! pftpar from LPJ-LMFire
! 1       2       3       4       5       6         7       8       9       10       11       12        13       14      15       16     17         18       19      20      21      22       23       24       25       26      27          28         29      30         31      32      33       34       35     36      37      38     39     40     41      42    43    44     45     46    47    48    49    50      51
! 0.850,  0.000,  0.000,  0.500,  0.100,  100.000,  2.000,  2.000,  20.000,  2.000,  29.000,  330.000,  29.000,  1.000,  1.000,  1.000,  1000.000,  30.000,  4.000,  3.000,  0.000,  2.000,   25.000,  30.000,  55.000,  0.900,  15.500,     1000.000,  0.000,  1000.000,  0.000,  0.333,  0.0301,  0.0281,  1.75,  -3.75,  2.52,  -0.78,  0.15,  0.12,  0.160,  0.5,  3.0,  0.05,  1580,  103,  6.8,  8.1,  8.5,  1.999,  15
! *0.700,  0.000,  0.350,  0.500,  0.100,  100.000,  0.500,  1.000,  20.000,  1.000,  29.000,  330.000,  29.000,  1.000,  3.000,  1.000,  1000.000,  30.000,  4.000,  3.000,  0.000,  2.000,   25.000,  30.000,  55.000,  0.900,  15.500,     1000.000,  0.000,  1000.000,  0.000,  0.100,  0.1085,  0.2120,  1.90,  -4.20,  3.40,  -2.11,  0.15,  0.50,  0.351,  0.5,  3.0,  0.40,  1664,  63,  2.2,  3.4,  8.5,  2.540,  15
! 0.700,  0.000,  0.000,  0.300,  1.000,  100.000,  2.000,  2.000,  20.000,  2.000,  29.000,  330.000,  29.000,  2.000,  1.000,  1.000,  1000.000,  30.000,  4.000,  3.000,  0.000,  -4.000,  20.000,  30.000,  42.000,  0.900,  -2.000,     22.000,  900.000,  1000.000,  0.000,  0.333,  0.0670,  0.5590,  2.57,  -6.20,  4.60,  -3.90,  0.15,  0.12,  0.094,  0.5,  3.0,  0.10,  1568,  106,  4.8,  5.7,  17.6,  3.240,  15
! 0.700,  0.000,  0.000,  0.500,  1.000,  100.000,  1.000,  1.000,  20.000,  1.000,  29.000,  330.000,  29.000,  1.000,  1.000,  1.000,  1000.000,  30.000,  4.000,  3.000,  0.000,  -4.000,  20.000,  30.000,  42.000,  0.800,  3.000,      18.800,  1200.000,  1000.000,  0.000,  0.333,  0.0451,  0.1412,  1.23,  -2.20,  6.86,  -7.30,  0.15,  0.50,  0.070,  0.5,  3.0,  0.10,  1568,  106,  4.8,  5.7,  17.6,  3.240,  15
! 0.800,  0.000,  0.000,  0.500,  1.000,  120.000,  0.500,  1.000,  20.000,  1.000,  29.000,  330.000,  29.000,  1.000,  2.000,  1.000,  200.000,   30.000,  4.000,  3.000,  0.000,  -4.000,  20.000,  25.000,  38.000,  0.800,  -17.000,    15.500,  1200.000,  1000.000,  0.000,  0.333,  0.0347,  0.1086,  1.68,  -3.53,  3.40,  -2.11,  0.15,  0.12,  0.094,  0.5,  3.0,  0.50,  1568,  106,  4.8,  5.7,  17.6,  3.240,  15
! *0.900,  0.000,  0.000,  0.300,  1.200,  100.000,  2.000,  2.000,  20.000,  2.000,  29.000,  330.000,  29.000,  2.000,  1.000,  1.000,  1000.000,  30.000,  4.000,  3.000,  1.000,  -4.000,  15.000,  25.000,  38.000,  0.800,  -32.500,    -2.000,  600.000,  23.000,  0.000,  0.333,  0.0292,  0.1086,  1.23,  -2.20,  7.95,  -8.92,  0.15,  0.12,  0.094,  0.5,  3.0,  0.44,  1568,  106,  4.8,  5.7,  17.6,  3.240,  15
! 0.900,  0.000,  0.000,  0.300,  1.200,  100.000,  0.500,  1.000,  20.000,  1.000,  29.000,  330.000,  29.000,  2.000,  2.000,  1.000,  200.000,   30.000,  4.000,  3.000,  1.000,  -4.000,  15.000,  25.000,  38.000,  0.900,  -1000.000,  -2.000,  350.000,  23.000,  0.000,  0.333,  0.0347,  0.1086,  1.32,  -2.46,  3.35,  -2.03,  0.15,  0.12,  0.094,  0.5,  3.0,  0.44,  1568,  106,  4.8,  5.7,  17.6,  3.240,  15
! 0.900,  0.000,  0.350,  0.500,  0.700,  100.000,  1.000,  1.000,  1.000,  2.000,   29.000,  0.000,    29.000,  3.000,  4.000,  0.750,  100.000,   0.000,   0.001,  3.000,  1.000,  -4.000,  10.000,  30.000,  45.000,  0.650,  -1000.000,  15.500,  0.000,  1000.000,  0.000,  -99,  -99,  -99,  -99,  -99,  -99,  -99,  0.15,  1.00,  -99.000,  0.5,  3.0,  0.50,  1568,  106,  4.8,  5.7,  17.6,  3.240,  2
! 0.900,  1.000,  0.350,  0.500,  0.200,  100.000,  1.000,  1.000,  1.000,  2.000,   29.000,  0.000,    29.000,  3.000,  4.000,  0.750,  100.000,   0.000,   0.001,  3.000,  0.000,  6.000,   20.000,  45.000,  55.000,  0.400,  15.500,     1000.000,  0.000,  1000.000,  0.000,  -99,  -99,  -99,  -99,  -99,  -99,  -99,  0.15,  1.00,  -99.000,  0.5,  3.0,  0.50,  1664,  63,  2.2,  3.4,  8.5,  2.540,  2
