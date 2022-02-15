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

use pftparmod,    only : npft,pftpar,tree
use statevarsmod, only : ndyear,dayvars,soilvars,gppvars,vegvars,topovars,lprint,gprint,gridlon,gridlat

integer(i4), intent(in) :: year
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day

!-------------------------
! Parameters
real(sp), parameter :: k       = 0.0548
real(sp), parameter :: tc      = 1. / 56.02   ! Constant in Arrhenius equation (Sitch et al. 2003; Eq. 23)
real(sp), parameter :: min_npp = 0.02

!-------------------------
! Pointer variables
real(sp), pointer :: tmean        ! 24 hour mean temperature (degC)
real(sp), pointer :: tsoil        ! Soil temperature (top layer) (K)
logical,  pointer, dimension(:) :: present
real(sp), pointer, dimension(:) :: gpp          ! Gross primary productivity under actual condition (g C m-2 d-1)
real(sp), pointer, dimension(:) :: npp0         ! Net primary productivity (g C m-2 d-1)
real(sp), pointer, dimension(:) :: npp_tot      ! Total net primary productivity (g C d-1)
real(sp), pointer, dimension(:) :: aresp        ! Autotrophic maintenence respiration (g C m-2 d-1)
real(sp), pointer :: dphen        ! Phenology status of summergreen (proportion of leaf-on) (fraction)
real(sp), pointer, dimension(:) :: lm_ind       ! Leaf carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: rm_ind       ! Root carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: sm_ind       ! Sapwood carbon mass of individual (gC m-2)
real(sp), pointer, dimension(:) :: nind         ! PFT population
real(sp), pointer, dimension(:) :: cstore       ! Carbon storage for growth respiration (gC m-2)
real(sp), pointer :: cellarea     ! Area of gridcell (m2)
real(sp), pointer :: areafrac     ! Ground area fraction in gridcell (fraction)

real(sp) :: tsoil_C               ! Soil temperature (degC)
real(sp) :: gtemp_air             ! Arrhenius equation varaible for air
real(sp) :: gtemp_soil            ! Arrhenius equation varaible for soil
real(sp) :: lresp_ind             ! Leaf respiration for individual (gC m-2 d-1)
real(sp) :: rresp_ind             ! Root respiration for individual (gC m-2 d-1)
real(sp) :: sresp_ind             ! Sapwood respiration for individual (gC m-2 d-1)
real(sp) :: aresp_ind             ! Total respiration for individual (gC m-2 d-1)

real(sp) :: respcoeff
real(sp) :: l_c2n
real(sp) :: r_c2n
real(sp) :: s_c2n

integer :: pft

!-------------------------

cellarea => topovars(grid)%cellarea
areafrac => topovars(grid)%areafrac

tmean   => dayvars(grid,day)%tmean

tsoil   => soilvars(grid)%Tsoil(2)

gpp     => gppvars(grid,day)%gpp
npp0    => gppvars(grid,day)%npp
npp_tot => gppvars(grid,day)%npp_tot
aresp   => gppvars(grid,day)%aresp

present => vegvars(grid)%present
dphen   => vegvars(grid)%dphen(day)
lm_ind  => vegvars(grid)%lm_ind
sm_ind  => vegvars(grid)%sm_ind
rm_ind  => vegvars(grid)%rm_ind
nind    => vegvars(grid)%nind
cstore  => vegvars(grid)%cstore

!-------------------------

if (areafrac <= 0.) then
  npp0 = -9999.
  npp_tot = 0.
  return
end if

!-------------------------

tsoil_C = tsoil - Tfreeze

tsoil_C = tmean ! Temporary fix for numerical instablity

! Calculate variables for Arrhenius equation (Sitch et al., 2003; Eq. 23)
if (tmean > -40.) then
  gtemp_air = exp(308.56 * (tc - 1. / (tmean  + 46.02)))
else
  gtemp_air = 0.
end if

if (tsoil_C > -40.) then
  gtemp_soil = exp(308.56 * (tc - 1. / (tsoil_C  + 46.02)))
else
  gtemp_soil = 0.
end if

!-------------------------

do pft = 1, npft

  if (present(pft)) then

    !------

    respcoeff = pftpar(5,pft)
    l_c2n     = pftpar(11,pft)
    r_c2n     = pftpar(13,pft)

    if (tree(pft)) then
      s_c2n = 1. / pftpar(12,pft)   ! To avoid divide by zero in grass PFTs
    else
      s_c2n = 0.
    end if

    !------

    ! Calculate respiratoin rate for each carbon compartment

    lresp_ind = (lm_ind(pft) / l_c2n) * respcoeff * k * gtemp_air  * dphen
    rresp_ind = (rm_ind(pft) / r_c2n) * respcoeff * k * gtemp_soil * dphen
    sresp_ind = (sm_ind(pft) * s_c2n) * respcoeff * k * gtemp_air

    ! Sum all maintenence respiration rate of individual
    aresp_ind = lresp_ind + rresp_ind + sresp_ind

    ! Multiply individual respiration with population to find total repiration of the day
    aresp(pft) = aresp_ind * nind(pft)

    if (isnan(aresp(pft))) write(0,*) 'NaN respiration at', gridlon(grid), gridlat(grid), nind,tmean,tsoil_C

    ! Subtract respiration from GPP to get NPP (gC m-2 d-1)
    npp0(pft) = gpp(pft) - aresp(pft)

    ! Make up negative NPP with carbon storage for growht respirtaion
    if (npp0(pft) < 0.) then
      npp0(pft) = npp0(pft) + cstore(pft)
      cstore(pft) = 0.
    end if

    ! Allocate 25% of npp to cstore, while reducing current cstore by 50%
    cstore(pft) = 0.25 * npp0(pft) + 0.5 * cstore(pft)

    cstore(pft) = max(cstore(pft), 0.)

    npp0(pft) = 0.75 * npp0(pft)

    ! Calculate total daily NPP of entire gridcell (gC d-1)
    npp_tot(pft) = npp0(pft) * cellarea * areafrac

    !-------------------------

    if (areafrac <= 0.) then
      npp0(pft) = -9999
      npp_tot(pft) = 0.
    end if

  end if

  ! if (lprint .and. grid==gprint .and. pft == npft) write(0,*) 'pft: ', year, day, npp0, vegvars(grid)%fpc_grid

end do ! PFT loop

!-------------------------

end subroutine npp

!---------------------------------------------------------------------

end module nppmod
