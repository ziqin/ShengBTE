# Files included in this package #

```
├── Doc
│   ├── arch.make.example
│   ├── ShengBTE.pdf
│   └── ShengBTE.tex
├── FILE_MAP
├── LICENSE
├── Src
│   ├── arch.make
│   ├── conductivity.f90
│   ├── config.f90
│   ├── data.f90
│   ├── input.f90
│   ├── integrals.f90
│   ├── iterations.f90
│   ├── Makefile
│   ├── misc.f90
│   ├── phonon_routines.f90
│   ├── processes.f90
│   ├── scaling.f90
│   ├── ShengBTE.f90
│   ├── symmetry.f90
│   └── wedgetc.f90
└── Test
    ├── CONTROL
    ├── FORCE_CONSTANTS_2ND
    ├── FORCE_CONSTANTS_3RD
    └── Reference
        ├── BTE.cumulative_kappa
        ├── BTE.cumulative_kappa_scalar
        ├── BTE.cumulative_kappa_tensor
        ├── BTE.cv
        ├── BTE.dos
        ├── BTE.kappa
        ├── BTE.kappa_scalar
        ├── BTE.kappa_sg
        ├── BTE.kappa_tensor
        ├── BTE.omega
        ├── BTE.P3
        ├── BTE.P3_total
        ├── BTE.pdos
        ├── BTE.qpoints
        ├── BTE.v
        ├── BTE.w
        ├── BTE.w_anharmonic
        ├── BTE.w_final
        └── BTE.w_isotopic
```

# Organization of the code #

The whole of ShengBTE is written in Fortran and contained in the Src subdirectory. What follows is a brief summary of the contents of each file:

* ShengBTE.f90: Main program

* data.f90: Physical constants, table of isotopic masses, and other relevant data

* misc.f90: Miscelaneous mathematical routines: cross product, integer division and Frobenius norm of a 3x3 matrix

* symmetry.f90: Thin wrapper around those parts of Atsushi Togo's [spglib](http://spglib.sourceforge.net/) relevant for ShengBTE. That library is used to detect the symmetry group of the system and obtain a matrix representation of its operations

* config.f90: Code to read and parse the CONTROL file and perform some operations based on this data. More specifically, this file contains calls to the symmetry routines to fill the associated data structures. A function used to compute the optimal broadening for each mode is also defined here on account of relying on information from CONTROL

* wedgetc.f90: The wedge() subroutine defined here harnesses the symmetries to reduce the gamma-centered regular grid used for sampling the phonon modes to a set of equivalence classes

* input.f90: Subroutines used to read FORCE\_CONSTANTS\_2ND and FORCE\_CONSTANTS\_3RD

* phonon_routines.f90: Code used to compute the phonon spectrum based on second-order IFC matrices and set of dielectric parameters

* integrals.f90: Routines used to calculate the lattice specific heat and the small-grain-limit reduced thermal conductivity

* processes.f90: Calculation of the number of allowed processed, their scattering rate and their phase space volume

* iterations.f90: Routines to set up the initial conditions for the iterative process and to perform iterations

* conductivity.f90: Code used to compute the thermal conductivity and its cumulative version once all the ingredients are known

* scaling.f90: Scaling of tau for nanowires oriented along an arbitrary direction of the bulk


# Other files #

* ShengBTE.tex, ShengBTE.pdf: User documentation for ShengBTE (LaTeX source and compiled version)

* LICENSE: Copy of the GNU General Public License, version 3

* Makefile: set of make rules for building ShengBTE

* arch.make.example: Machine-specific variables needed for the compilation. Must be copied to arch.make and customized appropriately

* CONTROL, FORCE\_CONSTANTS\_2ND, FORCE\_CONSTANTS\_3RD: Complete set of input files to run ShengBTE for InAs, a test case

* Reference/\*: Results obtained by the authors after running ShengBTE with the provided input files. Minor differences in the results can be expected depending on architecture, compiler and libraries

* FILE\_MAP.md: This file
