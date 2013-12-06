GLU
===

Gauge Link Utility (GLU) is a free, thread-parallel, number of dimensions (ND) and number of colors (NC) generic, lattice field
theory library for various gluonic observables. It is written in the c programming language.

The library can be used for:

File conversion : Support for ILDG, Scidac, NERSC, HiRep and MILC configuration files. Both reading and writing.

Gauge fixing : Landau and Coulomb gauge fixing is possible. With binding to the FFTW library, Fourier Accelerated routines
               for very fast gauge fixing have been implemented.

Smearing : APE, STOUT and LOG with their hypercubically blocked variants HYP, HEX and HYL smearing routines are available.
           So is the computation of the STOUT and LOG Wilson flow. Inclusion of various rectangle coefficients is
           available for the APE, STOUT and LOG and the Wilson flow routines.
           
Quenched U1 : The code allows for the computation of a SU(NC)xU(1) gauge field, where the U(1) gauge field is quenched.

Gauge correlators : Static potential computations, topological charge correlation functions and configuration space gluon 
                    propagators can be computed. With binding to FFTW, momentum space gluon correlation function 
                    calculations are also possible.
                    
Plus plenty more to be found in the documentation.
