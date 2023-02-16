module mpivarsmod

! Derived type containing job information to be sent to all processes

use parametersmod, only : i4,sp

implicit none

integer :: ofid

!-------------------------------------------------------

type mpivars
  character(200) :: cfile_spinup      ! infile name
  character(200) :: cfile_transient   ! infile name
  character(200) :: soilfile          ! infile name
  character(200) :: topofile          ! infile name
  character(200) :: outfile           ! outfile name
  character(200) :: gridlistfile      ! gridlist file name
  character(200) :: timestring        ! timestring with startyr/calcyr
  logical        :: dospinup
  logical        :: dotransient
  integer(i4)    :: spin_baseyr
  integer(i4)    :: spin_startyr
  integer(i4)    :: spinupyears
  integer(i4)    :: tran_baseyr
  integer(i4)    :: tran_startyr
  integer(i4)    :: transientyears
  integer(i4)    :: ilen              ! index length of validcells
  integer(i4)    :: spin_tlen
  integer(i4)    :: tran_tlen
  integer(i4)    :: outmode           ! output mode: 0 = global ncfile, 1 = grid textfile
  integer(i4)    :: nproc             ! number of processors

end type mpivars

!-------------------------------------------------------

end module mpivarsmod
