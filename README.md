GLU
===

Gauge Link Utility (GLU) is a free, thread-parallel, number of dimensions (ND) and number of colors (NC) generic, lattice field
theory library for various gluonic observables. It is written in the c programming language.

The library can be used for:

File conversion : Support for ILDG, Scidac, NERSC, HiRep, MILC, and CERN configuration files. Both reading and writing.

Gauge fixing : Landau and Coulomb gauge fixing is possible. With binding to the FFTW library, Fourier Accelerated routines
               for very fast gauge fixing have been implemented with the NLCG (default) or FASD algorithms.

Smearing : APE, STOUT and LOG with their hypercubically blocked variants HYP, HEX and HYL smearing routines are available.
           So is the computation of the STOUT and LOG Wilson flow. Inclusion of various rectangle coefficients is
           available for the APE, STOUT and LOG and the Wilson flow routines.
           
Quenched U1 : The code allows for the computation of a SU(NC)xU(1) gauge field, where the U(1) gauge field is quenched.

Gauge correlators : Static potential computations, topological charge correlation functions and configuration space gluon 
                    propagators can be computed. With binding to FFTW, momentum space gluon correlation function 
                    calculations are also possible.

Update : Heatbath updating and configuration generation is possible.
                    
Plus plenty more to be found in the documentation.

Compilation
===========

The usual,

    ./configure CC={} CFLAGS={} --prefix={}

If using the Fourier acceleration (I recommend you do) then

    --with-fftw={}

should be called.

If compiling for a different number of colors (NC) or dimensions (ND) (default is NC=3 and ND=4) the options

    --with-NC={} --with-ND={}

are your friends.

There a bunch of other options that are briefly synopsised at the top of configure.ac.

After that,

    make && make all install

will create the binary GLU in $prefix/bin/ and the library libGLU.a in $prefix/lib/ .

    make documentation

will create pdfs of the documentation in the $prefix/docs/ directory

    make doxygen

will create the doxygen with perhaps the callgraphs if the library dot is available.

the code either uses SSE2 intrinsics or defaults to generic code. I strongly recommend using the intrinsics branch (configure will tell you if it has found <immintrin.h> and openMP, this needs to be specified via the CFLAGS during compilation, so for example with icc:

    CFLAGS="-O3 -xHOST -fopenmp"
    
or with gcc

    CFLAGS="-O3 -msse2 -fopenmp"

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

Input File
==========

The input file expects a very specific format and even though you might not be using the gluon props, smearing, gauge-fixing, heatbath, whatever ... they still need to be there or the very simple reader will complain. Don't worry though, if they aren't used in your specific calculation they will just be ignored. 
