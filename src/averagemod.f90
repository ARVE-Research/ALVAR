module averagemod

implicit none

contains

subroutine longtermave(ayr,grid,day,ndyear)

use parametersmod, only : i4,sp,dp
use statevarsmod,  only : sv,ns,nl
use pftparmod,     only : npft

real(sp), intent(in) :: ayr
integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day
integer(i4), intent(in) :: ndyear

sv(grid)%vegvars%meanfpc = sv(grid)%vegvars%meanfpc + sv(grid)%vegvars%fpc_grid / ndyear / ayr

sv(grid)%longave%daet   = sv(grid)%longave%daet + sv(grid)%dayvars%daet(day) / ayr
sv(grid)%longave%gpp    = sv(grid)%longave%gpp + (sv(grid)%gppvars%gpp(day,1) + sv(grid)%gppvars%gpp(day,2) + &
                          sv(grid)%gppvars%gpp(day,6) &
                         + sv(grid)%gppvars%gpp(day,3) + sv(grid)%gppvars%gpp(day,4) + sv(grid)%gppvars%gpp(day,5) &
                         + sv(grid)%gppvars%gpp(day,7) + sv(grid)%gppvars%gpp(day,8) + sv(grid)%gppvars%gpp(day,9)) / ayr
sv(grid)%longave%npp    = sv(grid)%longave%npp + (sv(grid)%gppvars%npp(day,1) + sv(grid)%gppvars%npp(day,2) + &
                            sv(grid)%gppvars%npp(day,6) &
                         + sv(grid)%gppvars%npp(day,3) + sv(grid)%gppvars%npp(day,4) + sv(grid)%gppvars%npp(day,5) &
                         + sv(grid)%gppvars%npp(day,7) + sv(grid)%gppvars%npp(day,8) + sv(grid)%gppvars%npp(day,9)) / ayr
sv(grid)%longave%lm_ind = sv(grid)%longave%lm_ind + sum(sv(grid)%vegvars%lm_ind * sv(grid)%vegvars%nind) / ndyear / ayr
sv(grid)%longave%rm_ind = sv(grid)%longave%rm_ind + sum(sv(grid)%vegvars%rm_ind * sv(grid)%vegvars%nind) / ndyear / ayr
sv(grid)%longave%sm_ind = sv(grid)%longave%sm_ind + sum(sv(grid)%vegvars%sm_ind * sv(grid)%vegvars%nind) / ndyear / ayr
sv(grid)%longave%hm_ind = sv(grid)%longave%hm_ind + sum(sv(grid)%vegvars%hm_ind * sv(grid)%vegvars%nind) / ndyear / ayr
sv(grid)%longave%height = sv(grid)%longave%height + maxval(sv(grid)%vegvars%height) / ndyear / ayr
sv(grid)%longave%soilw  = sv(grid)%longave%soilw + &
                          (sum(sv(grid)%soilvars%Tliq(1:6)) / sum(sv(grid)%soilvars%Tsat(1:6))) / ndyear / ayr
sv(grid)%longave%Tsoil  = sv(grid)%longave%Tsoil + sum(sv(grid)%soilvars%Tsoil(1:3)) / 3. / ndyear / ayr
sv(grid)%longave%zsno   = sv(grid)%longave%zsno + sv(grid)%soilvars%zsno / ndyear / ayr


! print *, sv(grid)%longave%Tsoil-273.15, sv(grid)%longave%zsno, sv(grid)%longave%npp, sv(grid)%longave%soilw

end subroutine longtermave



subroutine output_ave(grid,day,ndyear)

use parametersmod, only : i4,sp,dp,Tfreeze
use statevarsmod,  only : sv,ns,nl
use pftparmod,     only : npft

integer(i4), intent(in) :: grid
integer(i4), intent(in) :: day
integer(i4), intent(in) :: ndyear


if (day == 1) sv(grid)%outvars%livebiomass = 0.
sv(grid)%outvars%livebiomass = sv(grid)%outvars%livebiomass + sum(sv(grid)%vegvars%lm_ind * sv(grid)%vegvars%nind) / ndyear &
                                                            + sum(sv(grid)%vegvars%rm_ind * sv(grid)%vegvars%nind) / ndyear &
                                                            + sum(sv(grid)%vegvars%sm_ind * sv(grid)%vegvars%nind) / ndyear

if (day == 1) sv(grid)%outvars%AET = 0.
sv(grid)%outvars%AET = sv(grid)%outvars%AET + sv(grid)%dayvars%daet(day)

if (day == 1) sv(grid)%outvars%GPP = 0.
sv(grid)%outvars%GPP = sv(grid)%outvars%GPP + sum(sv(grid)%gppvars%gpp(day,1:npft))

if (day == 1) sv(grid)%outvars%NPP = 0.
sv(grid)%outvars%NPP = sv(grid)%outvars%NPP + sum(sv(grid)%gppvars%npp(day,1:npft))

if (day == 1) sv(grid)%outvars%treecover = 0.
sv(grid)%outvars%treecover = sv(grid)%outvars%treecover + sum(sv(grid)%vegvars%fpc_grid(1:7)) / ndyear

if (day == 1) sv(grid)%outvars%grasscover = 0.
sv(grid)%outvars%grasscover = sv(grid)%outvars%grasscover + sum(sv(grid)%vegvars%fpc_grid(8:9)) / ndyear

if (day == 1) sv(grid)%outvars%cover = 0.
sv(grid)%outvars%cover = sv(grid)%outvars%cover + sv(grid)%vegvars%fpc_grid(1:npft) / ndyear

if (day == 1) sv(grid)%outvars%nind = 0.
sv(grid)%outvars%nind = sv(grid)%outvars%nind + sv(grid)%vegvars%nind(1:npft) / ndyear

if (day == 1) sv(grid)%outvars%crownarea = 0.
sv(grid)%outvars%crownarea = sv(grid)%outvars%crownarea + sv(grid)%vegvars%crownarea(1:npft) / ndyear

if (day == 1) sv(grid)%outvars%soilmoisture = 0.
sv(grid)%outvars%soilmoisture = sv(grid)%outvars%soilmoisture + sv(grid)%soilvars%Tliq(1:nl) / sv(grid)%soilvars%Tsat(1:nl) / ndyear

if (day == 1) sv(grid)%outvars%soiltemp = 0.
sv(grid)%outvars%soiltemp = sv(grid)%outvars%soiltemp + (sv(grid)%soilvars%Tsoil(1:nl) - Tfreeze) / ndyear

if (day == 1) sv(grid)%outvars%height = 0.
sv(grid)%outvars%height = sv(grid)%outvars%height + sv(grid)%vegvars%height(1:npft) / ndyear

if (day == 1) sv(grid)%outvars%GDD5 = 0.
sv(grid)%outvars%GDD5 = sv(grid)%outvars%GDD5 + max((sv(grid)%dayvars%tmean(day) - 5.),0.)


end subroutine output_ave




end module averagemod
