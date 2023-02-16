module rootmod

! Module to calculate soil layer root fraction
! Code adapted from ARVE-DGVM by Leo Lai (Oct 2022)

implicit none

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine rootdist(tree,present,zipos,lm_ind,rm_ind,sm_ind,hm_ind,nind,crownarea,frac,rootdepth,rootfracl)

use parametersmod, only : i4,sp,dp
use pftparmod,     only : npft,pftpar
use statevarsmod,  only : nl,ns,soildepth

implicit none

logical,  dimension(:),   intent(in)    :: tree         ! Tree PFT
logical,  dimension(:),   intent(in)    :: present      ! PFT is present in grid cell
real(sp), dimension(ns:nl+1), intent(in)   :: zipos          ! Snow/soil layer interface z position (m)
real(sp), dimension(:),   intent(in)    :: lm_ind       ! Individual leaf mass (gC)
real(sp), dimension(:),   intent(in)    :: rm_ind       ! Individual fine root mass (gC)
real(sp), dimension(:),   intent(in)    :: sm_ind       ! Individual sapwood mass (gC)
real(sp), dimension(:),   intent(in)    :: hm_ind       ! Individual heartwood mass (gC)
real(sp), dimension(:),   intent(in)    :: nind         ! Gridcell individual density (indiv/m2)
real(sp), dimension(:),   intent(in)    :: crownarea    ! Tree crownarea per individual (m2)
real(sp), dimension(:),   intent(inout) :: frac         ! Fraction of sapwood and heartwood that is below ground (in coarse roots)(gC)
real(sp), dimension(:),   intent(inout) :: rootdepth    ! Root depth (m)
real(sp), dimension(:,:), intent(inout) :: rootfracl    ! Root fraction in that soil layer
! real(dp), pointer, dimension(:,:) :: rm_max           ! Maximal individual fine root mass (gC)

real(sp), parameter :: frac_cr_tot = 0.2        ! Percent of the above ground biomass that goes to coarse roots (see refs in Ryan et al.
real(sp), parameter :: OMorgC      = 0.5800464  ! g org C =sorg/ 1.724 (conversion factor from Nelson & Sommers 1996)

!-------------------------
! Local variables
real(sp) :: abar                                ! Parameter representing mean root distribution profile
real(sp) :: Bbar                                ! Average standing root biomass (kg m-2)
real(sp) :: littleb                             ! Parameter representing variable root distribution profile (=abar * Bbar exp alpha)
real(sp) :: Broot                               ! Root biomass in kg organic/m2
real(sp) :: rootgro                             ! Parameter for root growth direction (0 -1 )
real(sp) :: Bunit
real(sp) :: a                                   ! Vegetation dependent inverse e-folding length scale
real(sp) :: b
real(sp) :: trm_ind                             ! Individual total root mass (gC)
real(sp) :: term_a                              ! Temporary variable
real(sp) :: term_b                              ! Temporary variable
real(sp), dimension(nl) :: rootf                ! Individual cumulative distribution profiles
!real(dp), dimension(nl,npft) :: rootdens       ! Root density per pft and per layer (g/m2)
integer(i4) :: l
integer(i4) :: pft

!-----------------------------------------------------------

do pft = 1, npft

  if (present(pft)) then

    abar    = pftpar%abar(pft)
    Bbar    = pftpar%Bbar(pft)
    littleb = pftpar%littleb(pft)

    !-------------------------
    ! The fraction of sapw and heartw below ground is based on the assumption that coarse roots
    ! are 20% of total aboveground biomass (see Ryan, Binkley & Formes Adv.Ecol. Res. 1997)
    if (tree(pft)) then

      frac(pft) = (0.2 * (lm_ind(pft) + hm_ind(pft) + sm_ind(pft))) / (hm_ind(pft) + sm_ind(pft))
      trm_ind   = rm_ind(pft) + frac(pft) * (sm_ind(pft) + hm_ind(pft))

    else

      trm_ind = rm_ind(pft)

    end if

    !-------------------------
    ! Find the rooting depth

    ! Intialize root growth direction (0.80 value is from Arora & Boer)
    ! If we wish to add in an ability for the allocation of resources to change dependent upon
    ! something like PET, climate, or something else, we would change rootgro here. It is the
    ! determinant of whether the roots go preferentially vertical or horizontally. 1 = vertical, 0 = horizontal
    ! if try to put in climate, grasses are heavily swayed by mean annual precip, trees are not. See
    ! Schenk & Jackson J Ecol. 2002 v90 480-494. However with the daily C allocation, this is effectively occuring.
    rootgro = 0.8

    ! Put root biomass (gC/indiv) into kg organic/m2 of vegetated surface
    Broot = trm_ind * 1.e-3  / crownarea(pft) / OMorgC

    ! Find the rooting depth
    rootdepth(pft) = 3. * (Broot**rootgro) / littleb

    if (rootdepth(pft) > soildepth) then

      ! If the rooting depth is deeper than the soil depth the root distribution profile is not allowed to deepen
      ! with increasing root biomass. This is equivalent to reducing rootgro since to keep rootdepth constant at soil
      ! depth it requires that soildepth = 3/abar (Broot/Bbar)exp (rootgro) hence rootgro =
      ! ln(soildepth*abar/3)/ln(Broot/Bbar)

      ! NOTE: the parameterization from the paper does not make much sense. It says that the rootgro value should
      ! decrease as the roots get to the max soil depth but the value actually increases, exponentially. It is easy
      ! to plot and check. Instead of that relation I am just setting alpha to 0 when the rootdepth is greater than
      ! the soil depth. This makes all future growth be horizontal as you would expect (this assumes a constant soil
      ! depth across a grid cell and no cracks in the bedrock for roots to exploit

        ! rootgro  = log(soildepth * abar(pft) / 3._dp) / log(Broot / Bbar(pft)) Arora pub relation
        rootgro = 0.

        ! Now set the rooting depth to be the soil depth
        rootdepth(pft) = soildepth

    end if

    !-------------------------

    b = abar * (Bbar**rootgro)

    a = b / (Broot**rootgro)

    do l = 1, nl

     ! These calculations can produce underflow errors, as they exponentially decline.
     ! this is limited by these max statements below. JM 10.11.2010
     ! Update: This was likely a result of the problem with root depth once below soil depth, so likely not needed anymore. JM 17.03.2011

     ! Calculate the individual cumulative distribution profiles  (eqn 16)
     term_a   = max(-50.,(-abar * (Bbar / Broot)**rootgro * zipos(l+1)))
     rootf(l) = 1. - exp(term_a)

     ! Find the root density at this soil depth  !Eqn 9
     term_b = max(-50., (-a * zipos(l+1)))
     Bunit  = exp(term_b)

     ! Root density, can be useful for diagnostics but not strictly needed.
     ! rootdens(l,pft) = b * Broot**(1._dp - rootgro) * Bunit !* 1000._dp !convert to g/m2

   end do ! Soil layers loop

   !-------------------------
   ! Make rootfraction per layer, not cumulative
   rootfracl(pft,1) = rootf(1)

    do l = 2, nl
      rootfracl(pft,l) = rootf(l) - rootf(l-1)
    end do

    if (rootfracl(pft,1) == 1._dp) then  !don't allow all of root mass to be in the upper layer. Spread it between the top 2
        rootfracl(pft,2) = 0.5
        rootfracl(pft,2) = 0.5
    end if

    !-------------------------

  end if ! IF present loop

end do ! PFT loop


end subroutine rootdist


end module rootmod
