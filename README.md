Version 1.0 of ALVAR model (Leo O Lai, HKU, 28 Apr 2021)
Code adapted from gwgen weather generation, LPJ_LMFire and ARVE-DGVM (from ARVE-Research Git Repositories)

! Before compiling on epscc
$ cd ALVAR/
$ module load netCDF-Fortran/4.5.3-gmpich-2021.01 Autotools pkgconfig

! To compile:
$ ./autogen.sh
$ ./configure
$ make

! To run: program was written to run ALL validcells globally, but number of years can be specified
$ module purge              ! because Autotools switched GCC back to 9.3.0. v10.2 is required to run the code
$ module load netCDF-Fortran/4.5.3-gmpich-2021.01

$ mpirun -np 18 ./src/alvar /home/terraces/datasets/dgvm_input/climate/transient1871-2010_list-formatted.nc -110.25/40.25 1990/10 test.nc

! NOTE: user-specified lon/lat must be **EXACT of the input file values** (routine to find closest gridcell is yet to be written) - Leo Lai 28 Apr 2021)
