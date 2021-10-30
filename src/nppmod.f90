module nppmod

! Module to calculate NPP coded for ALVAR
! Code adapted from LPJ-LMFire by Leo Lai (Oct 2021)

use parametersmod, only : i2,i4,sp,dp,missing_i2,missing_sp,Tfreeze,daysec

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
real(sp), parameter :: respcoeff = 0.66
real(sp), parameter :: tc = 1. / 56.02
real(sp), parameter :: l_c2n = 29
real(sp), parameter :: r_c2n = 29
real(sp), parameter :: s_c2n = 330

!-------------------------
! Pointer variables
real(sp), pointer :: tmean
real(sp), pointer :: tsoil
real(sp), pointer :: lm_ind
real(sp), pointer :: rm_ind
real(sp), pointer :: sm_ind
real(sp), pointer :: nind
real(sp), pointer :: gpp0
real(sp), pointer :: npp0
real(sp), pointer :: npp_tot
real(sp), pointer :: aresp
real(sp), pointer :: cellarea
real(sp), pointer :: areafrac

real(sp) :: gtemp_air
real(sp) :: gtemp_soil
real(sp) :: lresp_ind
real(sp) :: rresp_ind
real(sp) :: sresp_ind
real(sp) :: aresp_ind

!-------------------------

tmean => dayvars(grid,day)%tmean
tsoil => soilvars(grid)%Tsoil(1)
gpp0 => vegvars(grid,day)%gpp
npp0 => vegvars(grid,day)%npp0
npp_tot => vegvars(grid,day)%npp_tot
aresp => vegvars(grid,day)%aresp
lm_ind => vegvars(grid,day)%lm_ind
sm_ind => vegvars(grid,day)%sm_ind
rm_ind => vegvars(grid,day)%rm_ind
nind => vegvars(grid,day)%nind

cellarea => topovars(grid)%cellarea
areafrac => topovars(grid)%areafrac


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


lresp_ind = respcoeff * (lm_ind / l_c2n) * gtemp_air  !* dphen(:,i)
rresp_ind = respcoeff * (rm_ind / r_c2n) * gtemp_soil !* dphen(:,i)
sresp_ind = respcoeff * (sm_ind / s_c2n) * gtemp_air

aresp_ind = lresp_ind + rresp_ind + sresp_ind

aresp = aresp_ind * nind

npp0 = gpp0 - aresp

npp_tot = npp0 * cellarea * areafrac

if (areafrac == 0.) then
  npp0 = -9999
  npp_tot = 0.
end if

! if (aresp > 0. .and. npp0 > 0.) print *, gpp0, aresp, npp0

end subroutine npp

!---------------------------------------------------------------------

end module nppmod
