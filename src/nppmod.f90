module nppmod

! Module to calculate NPP coded for ALVAR
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine npp(year,grid,day)

use statevarsmod, only : ndyear,dayvars,soilvars,vegvars,topovars

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: respcoeff = 0.66         ! Respiration coefficient (gC gN-1 d-1)
real(sp), parameter :: tc        = 1. / 56.02   ! Constant in Arrhenius equation (Sitch et al. 2003; Eq. 23)
real(sp), parameter :: l_c2n     = 29.          ! Leaf C:N ratio (Sitch et al. 2003; Table 3)
real(sp), parameter :: r_c2n     = 29.          ! Root C:N ratio (Sitch et al. 2003; Table 3)
real(sp), parameter :: s_c2n     = 330.         ! Sapwood C:N ratio (Sitch et al. 2003; Table 3)

!-------------------------
! Pointer variables
real(sp), pointer :: tmean        ! 24 hour mean temperature (degC)
real(sp), pointer :: tsoil        ! Soil temperature (top layer) (K)
real(sp), pointer :: gpp          ! Gross primary productivity under actual condition (g C m-2 d-1)
real(sp), pointer :: npp0         ! Net primary productivity (g C m-2 d-1)
real(sp), pointer :: npp_tot      ! Total net primary productivity (g C d-1)
real(sp), pointer :: aresp        ! Autotrophic maintenence respiration (g C m-2 d-1)
real(sp), pointer :: lm_ind       ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer :: rm_ind       ! Root carbon mass of individual (gC m-2)
real(sp), pointer :: sm_ind       ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer :: nind         ! PFT population
real(sp), pointer :: cellarea     ! Area of gridcell (m2)
real(sp), pointer :: areafrac     ! Ground area fraction in gridcell (fraction)

real(sp) :: tsoil_C               ! Soil temperature (degC)
real(sp) :: gtemp_air             ! Arrhenius equation varaible for air
real(sp) :: gtemp_soil            ! Arrhenius equation varaible for soil
real(sp) :: lresp_ind             ! Leaf respiration for individual (gC m-2 d-1)
real(sp) :: rresp_ind             ! Root respiration for individual (gC m-2 d-1)
real(sp) :: sresp_ind             ! Sapwood respiration for individual (gC m-2 d-1)
real(sp) :: aresp_ind             ! Total respiration for individual (gC m-2 d-1)

!-------------------------

tmean   => dayvars(grid,day)%tmean
tsoil   => soilvars(grid)%Tsoil(1)
gpp     => vegvars(grid,day)%gpp
npp0    => vegvars(grid,day)%npp
npp_tot => vegvars(grid,day)%npp_tot
aresp   => vegvars(grid,day)%aresp
lm_ind  => vegvars(grid,day)%lm_ind
sm_ind  => vegvars(grid,day)%sm_ind
rm_ind  => vegvars(grid,day)%rm_ind
nind    => vegvars(grid,day)%nind

cellarea => topovars(grid)%cellarea
areafrac => topovars(grid)%areafrac

!-------------------------

tsoil_C = tsoil - Tfreeze

! Calculate variables for Arrhenius equation (Sitch et al., 2003; Eq. 23)
if (tmean > -40.) then
  gtemp_air = exp(308.56 * (tc - 1. / (tmean  + 46.02)))
else
  gtemp_air = 0.
end if

if (tsoil > -40.) then
  gtemp_soil = exp(308.56 * (tc - 1. / (tsoil  + 46.02)))
else
  gtemp_soil = 0.
end if

!-------------------------
! Calculate respiratoin rate for each carbon compartment

lresp_ind = respcoeff * (lm_ind / l_c2n) * gtemp_air  !* dphen(:,i)
rresp_ind = respcoeff * (rm_ind / r_c2n) * gtemp_soil !* dphen(:,i)
sresp_ind = respcoeff * (sm_ind / s_c2n) * gtemp_air

! Sum all maintenence respiration rate of individual
aresp_ind = lresp_ind + rresp_ind + sresp_ind

! Multiply individual respiration with population to find total repiration of the day
aresp = aresp_ind * nind

! Subtract respiration from GPP to get NPP (gC m-2 d-1)
npp0 = gpp - aresp

! Calculate total daily NPP of entire gridcell (gC d-1)
npp_tot = npp0 * cellarea * areafrac

!-------------------------

if (areafrac == 0.) then
  npp0 = -9999
  npp_tot = 0.
end if

!-------------------------

! if (aresp > 0. .and. npp0 > 0.) print *, gpp, aresp, npp0

end subroutine npp

!---------------------------------------------------------------------

end module nppmod
