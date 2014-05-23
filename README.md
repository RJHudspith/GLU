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

Compilation
===========

The usual,

./configure --prefix={}

If using the Fourier acceleration then

--with-fftw={}

should be called.

If compiling for a different number of colors (NC) or dimensions (ND) the options

--with-NC={} --with-ND={}

are your friends.

If you want to include the thread-parallel FFTW routines then

CFLAGS="-fopenmp" --enable-OMP_FFTW

Should be used, this then looks in the FFTW directory that you have compiled for
libfftw3_omp.a

There a bunch of other options that are briefly synopsised at the top of configure.ac.

After that,

make && make all install

will create the binary GLU in $prefix/bin/ and the library libGLU.a in $prefix/lib/ .

make documentation

will create pdfs of the documentation in the $prefix/docs/ directory

make doxygen

will create the doxygen with perhaps the callgraphs if dot is available.

Usage
=====

Standard usage is

./GLU -i {input_file} -c {configuration file} -o {output file}

The output file is not necessary, the other two are. An example input file should be in $prefix/bin/
a couple of example input files can be printed to stdout by the command,

./GLU --autoin={COULOMB,LANDAU,STATIC_POTENTIAL,SUNCxU1,WFLOW}

If you would like more information about some option in the input file,

./GLU --help={option}

Should provide you with some information.
