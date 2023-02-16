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




end module averagemod
