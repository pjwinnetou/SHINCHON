MUSIC  (beta version)             {#mainpage}
===================

This is a beta version of the incoming public version of the heavy ion
hydronamics code MUSIC. 

MUSIC requires GSL libraries to be installed from
[here](https://www.gnu.org/software/gsl/)
as well as MPI (e.g. MPICH or OpenMPI) and a C++ compiler.  It only runs
properly on POSIX operating systems (tested on LINUX and Mac OS X).

---

## Compiling the code

* The MUSIC code can be compiled using standard cmake. 

Under the root directory, type the following commands:

mkdir build

cd build

cmake ..

make

make install

Once these commands are done, an executable "mpihydro" will be generated in the
root folder.

* Alternatively, you can compile the MUSIC code using a makefile.  In this
  option, please make sure the information about the MPI compiler and the
directory of the GSL library is correct in the src/GNUmakefile.  Then one can
compile the code by typing `make clean` and then `make` under the root folder.

## Run MUSIC with an input file

MUSIC is usually run under the root directory.

An input file is required that contains the line "EndOfData", preceded by a list
of parameter names and values, one per line, with parameter names and values
separated by a space or tab.  If omitted, each parameter will be assigned a
default value.  The list of possible parameters and default values, along with a
brief description, can be found in "utilities/music_parameters_dict.py" and in
the MUSIC documentation ("notes.pdf")

To run hydro simulation, one can start by running the python script
"utilities/generate_music_inputfile.py" to generate input files and job running
script. The job script can run under bash and qsub systems.  For some help
information about "utilities/generate_music_inputfile.py", one can simply type:
`generate_music_inputfile.py -h`

### Run MUSIC on multiple CPU cores

MUSIC can be invoked by MPI.  For example, to run on two processors and use the
sample input file, type: `mpirun -n 2 ./mpihydro input_example`

## Utilities

There are a few useful python and bash scripts stored in the utilities/ folder.
They are helpful for large-scale event-by-event numerical simulations.

---

## Publishing results computed with MUSIC

MUSIC is a 3+1D relativistic fluid dynamics simulation that can model the
evolution of the plasma formed in heavy ion collisions, compute the thermal
spectra of hadrons using Cooper-Frye, evaluate the decay of unstable hadrons,
and compute various hadronic observables with the thermal or post-decay spectra.
New results computed with MUSIC may necessitate modifying certain parts of the
code (routines computing hadronic observables, transport coefficients, initial
conditions, ...). To avoid any confusion about the physics behind MUSIC, we ask
that publications containing results computed with MUSIC specify the initial
conditions, equation of state, transport coefficients (first and second order),
freeze-out criteria and other information necessary to reproduce the results.

The following papers should be cited when referring to MUSIC:
1) B. Schenke, S. Jeon, C. Gale. "3+1D hydrodynamic simulation of relativistic heavy-ion collisions" Phys.Rev.C 82, 014903 (2010) [arXiv:1004.1408]
2) B. Schenke, S. Jeon, C. Gale. “Elliptic and triangular flow in event-by-event (3+1)D viscous hydrodynamics” Phys.Rev.Lett. 106, 042301 (2011) [arXiv:1009.3244]
3) J.-F. Paquet, C. Shen, G. S. Denicol, M. Luzum, B. Schenke, S. Jeon and C. Gale. "Production of photons in relativistic heavy-ion collisions" Phys. Rev. C 93, 044906 (2016) [arXiv:1509.06738]

## Modifying MUSIC

MUSIC users may need to modify various parts of the code. The GPL licence allows
for such modifications to be made and distributed (subject to restrictions that
can be found in the licence). Nevertheless, we ask and require that you state
clearly any scientifically relevant modifications made to the code when you
share this modified version of MUSIC. This is to avoid any misunderstanding
about the physics behind the code.

## Contributing to MUSIC

Should you make changes to MUSIC that you believe would benefit from being
shared with the wider community, please contact the code’s maintainers through
the MUSIC website. Contributions will be considered whenever possible. While no
promises can be made regarding the inclusion of these modifications in the main
version of the code, they could be made available separately on the MUSIC
website for other users’ convenience. Bug fixes will be considered separately
and will be applied promptly.
