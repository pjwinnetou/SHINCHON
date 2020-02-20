B3D v2.3

Adding higher order cumulants by Abbas Hassani
=================================
Important file for cumulants:

	  madai_repo/rhic/anal/src/anal_qn.cc

Look for 'Modifications' or 'AH,MH'
Copy 'anal_qn.cc' to  'anal_qn_ho.cc' and implement the h.o. cumulants; 
Change function names for 'anal_qn_ho.cc';
Output the 2n-order cumulants (already implemented) and the 4th-, 6th-order cumulants into different files;


Parallelisation by Ryan Duane Weller
====================================
Important file for parallelisation:

	  madai_repo/rhic/anal/src/anal_hbt.cc

Look for 'Starting parallelisation' or 'RDW,MH'
There are a few lines with pseudo-code; probably the OpenMP parallelisation is faster to do; however, if one knows which variables have to be passed on to the threads, then MPI parallelisation is also possible, but still more labour-intensive.


Known problems
==============
Running 'b3d' or 'analyze' requires this hack

	export HDF5_DISABLE_VERSION_CHECK=1

The executables report that the headers are version HDF5 1.8.10, but the library is 1.8.15 (current verison on Eridanus).
HDF5 is not part of the code package, so compiling should use the same headers and library version
Conjecture 1: maybe some files are pre-compiled and, thus, they are not updated when doing 'make'
Conjecture 2: HDF5 is not properly installed and there are still some residues of 1.8.10. Unlikely, because similar problem on a bare bone cluster with only local installs persists as well.


Compiling (legacy entry)
========================

* You need to modify the directory information in "madai_repo/makefile_defs.mk"
(e.g. path to madai_repo and hdf5 libraries)

* You need to update the symbolic link "madai" to point to "install/progdata/madai"

NOTE: you will need hdf5 libraries (these are not provided with the b3d.tgz file) on your system!!!!

Acutal compiling
================

Run

      make

to create the excecutables 'b3d' and 'analyze';

Clean up via
      
      make clean

Running the code	
================
Need to be run in order:

1. B3D routines (parameter after executable is important!) in the main directory
./b3d default

2. Analyze routines (parameter after executable is important!) in the main directory
./analyze default


Miscellaneous
=============

Added thesis by Ante Bilandzic and another paper on flow cumulants into the main directory.
