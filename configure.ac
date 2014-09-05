## GLU Automake instructions outline ##

# --prefix=  :: tell it where to install the binary
# --with-NC= :: number of colors, if not specified defaults to 3
# --with-ND= :: number of dimensions, if not specified defaults to 4
# --enable-single :: working precision is single, default is double
# --enable-OMP_FFTW :: enable OMP'd fftw routines (FFTW v3.3 and greater)
# --with-fftw= :: directory for the installed FFTW routines if --with-fftw
#                 argument is specified without a directory, will look for
#                 -lfftw3(f) in /usr/lib/ 
# --with-gsl= :: directory for GSL
# --with-lapacke= :: directory for lapacke install, used for large-NC eigen
#                    operations, good binding with complex.h
# --with-deriv= :: derivative type in gauge fixing options are :
#                  deriv_linn , deriv_full , deriv_fulln (last two are
#  		   exact logarithm variations)
# --enable-exp_exact :: exact exponentiation into the group being used for the
#                       gauge transformation matrices in the gauge fixing
# --enable-fastsmear :: dangerous, inaccurate but quick smearing
# --enable-notcondor :: enable the single-machine cache-ing and architecture
#                       dependent optimisations
# --with-GRACEY_PROJ= :: MOMggg scheme projection number, 0 is default
# --with-TOP_VALUE= :: Specify number of smearing steps before measuring the
#                      topological charge 
# --enable-GIVENS_APE :: use the Cabbibo-Marinari trace-maximisation for
#                        the APE and HYP smearings
# --enable-GLU_bgq :: inlined matrix multiplies via macros
# --with-IMPROVED_C0=  :: leading 1x1 term if not using IWASAKI, SYMANZIK or DBW2
# --with-IMPROVED_C1=  :: subleading 2x1 term only effective if C0 is defined
# --with-IMPROVED_STAPLE= :: Rectangle coefficients, choices are :
# 			     IWASAKI, SYMANZIK and DBW2
# --with-OVERIMPROVED_EPSILON= :: specify a value for the overimprovement
#                                 factor, default is 0.0
# --enable-verbiage :: more verbose output than usual
# --enable-ascii_check :: prints to stdout the data we would write to our
#                         binary files
# --enable-CLOVER_IMPROVE :: Are we using the O(a^4) clover def?
# --with-CLOVER_K5= :: specify the multiplicative term for the 3x3 wilson
#                      loop for the heavily improved field strength tensor
# --with-RNG= :: choices are KISS, MWC and GSL. Default is the WELL
# --with-LUXURY_GAUGE= :: choose the "best" (one that minimises the
# 		       	  functional) from specified number of random
#                         initial transformations, default is 20
# --enable-GLU_GFIX_SD :: turn on the steepest descent codes
# --enable-OVERRELAX_GF :: switch to the overrelaxed codes for gauge fixing

AC_INIT([GLU], [1.1], [renwick.james.hudspith@gmail.com])
AM_INIT_AUTOMAKE([-Wall])
AC_PROG_CC
AM_PROG_CC_C_O
AC_LANG([C]) 
AC_PROG_RANLIB ## include a check for libtool if not use ranlib?
AC_C_BIGENDIAN 

## SET up the m4
AC_CONFIG_MACRO_DIR([m4])

## I strongly recommend these flags for compiling with gcc
#CFLAGS="-O3 -mfpmath=sse -march=native -fopenmp"

## prefix!!
AC_DEFINE_UNQUOTED([HAVE_PREFIX],["$prefix"],[Look for where we prefix to])
	  
## Tells us where we are prefixed to
AC_MSG_NOTICE([Configuring for location "$prefix"])

## If we prefix away from the default /usr/ then create binary GLU
## will still try and create shared object libraries there
AM_CONDITIONAL([PREF],[ test "$prefix" == "NONE" ])
if test "$prefix" == "NONE" ; 
then 
    AC_MSG_NOTICE([Not creating binary "GLU"])
else
    AC_MSG_NOTICE([Creating binary "GLU" in $prefix/bin])
    ## check for doxygen
    AC_CHECK_PROG([DOXYGEN],[doxygen],["yes"],["NONE"])
    AM_CONDITIONAL([DOX],[ test "$DOXYGEN" == "NONE" ])
    ## check for dot
    AC_CHECK_PROG([GRAPHVIZ],[dot],["yes"],["NONE"])
    AM_CONDITIONAL([DOT],[ test "$GRAPHVIZ" == "NONE" ])
    ## check for pdflatex
    AC_CHECK_PROG([PDFLATEX],[pdflatex],["yes"],["NONE"])
    AM_CONDITIONAL([PDFLTX],[ test "$PDFLATEX" == "NONE" ])
    mkdir $prefix  
    mkdir $prefix/bin
    mkdir $prefix/docs/ ## new documentation folder ...
    ## copy the input file and a couple of small tests to bin	
    cp -v ./CONF/* $prefix/bin/
fi
AC_SUBST(PREF)

## Compile for a set NC, at the moment NC=2,3 supported may change 
AC_ARG_WITH([NC],
	AC_HELP_STRING([--with-NC=],
		       [Compile in the gauge group number NC]),
            [
	    AC_MSG_NOTICE([User specified NC, compiling for SU(${with_NC})])
	    AC_DEFINE_UNQUOTED([NC], [${with_NC}] , [Compiled for SU(NC)] )	
	    ],[
	    AC_MSG_NOTICE([User unspecified NC, default to SU(3)])	
	    ])

## C standard libraries
AC_CHECK_HEADERS([math.h],[
LDFLAGS="$LDFLAGS -lm" ])

## C standard library otherwise (dodgily) do it ourself
AC_CHECK_HEADERS([stdint.h],[])

## check for the unix memory stuff,
AC_CHECK_HEADERS([unistd.h], 
	[AC_MSG_NOTICE([We are attempting to be memory-safe in our operation.])],[])

## Check for sys/time.h
AC_CHECK_HEADERS([sys/time.h],
	[AC_MSG_NOTICE([We are intrinsically timing our operations.])],[])

## Check for the time.h, which we use for the date
AC_CHECK_HEADERS([time.h],
	[AC_MSG_NOTICE([We are providing an accurate date.])],[])

## I doubt we will not have these as they are all c-standard, check anyway
AC_CHECK_HEADERS([complex.h strings.h stdio.h stdlib.h limits.h],[],[
AC_MSG_ERROR([Sorry, but we need these headers (complex.h , strings.h , stdlib.h , stdio.h and limits.h) somewhere in the code])])

## And we can now run in single precision
single=false
AC_ARG_ENABLE([single],
	[ --enable-single Have our links in single or double precision],
	[
	case "${enableval}" in
	yes) single=true ## needed for fftw
	     AC_MSG_NOTICE([Using single precision storage for gauge fields]) 
	     AC_DEFINE([SINGLE_PREC],[],[Storage of single precision fields and utilisation of fftwf routines.]) 
	     ;;
	no) AC_MSG_NOTICE([Using double precision storage for gauge fields (default)]) ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-single]) ;;
	esac
	],[])

######################################### FFTW3 ###############################################
#                                                                                             #
#    This one is tricky, looks first for the static link to FFTW library, if not found        #
#    then the shared. If --with-fftw= then it looks in /usr/lib for the shared                #
#    otherwise we miss out on all of the Fourier transform stuff like Fourier acceleration,   #
#    Momentum space gluons and quenched U(1). Not a necessity.                                #
#                                                                                             #
###############################################################################################

## I have now written the code to utilise OpenMP parallel FFT routines
AC_ARG_ENABLE([OMP_FFTW],
	[ --enable-OMP_FFTW Utilise FFTW's OpenMP parallel FFTs],
	[case "${enableval}" in
	yes) OMP_FFTW=true ;;
	no) OMP_FFTW=false ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-OMP_FFTW]) ;;
	esac],[OMP_FFTW=false])

## Look for openmp, cannot turn on the omp fftw routines without it!
AC_CHECK_HEADERS([omp.h],[],[OMP_FFTW=false])

## Look for fftw in the useer specified place
fftw=false 
AC_ARG_WITH([fftw],
	AC_HELP_STRING([--with-fftw=<name>],
		       [Specify FFTW location for the GF]),
            [
	    ## see if we already have it	
	    AC_CHECK_HEADERS([fftw3.h],[local_fftw=true],[local_fftw=false])
	    
	    FFTW=${with_fftw}

	    ## are we compiling for single precision or not?
	    if test x$single = xtrue ; then
	       float_version=f
	    else	
	       float_version=""
	    fi

	    AC_CHECK_HEADERS([${FFTW}include/fftw3.h],[
	    AC_MSG_NOTICE([User defined fftw library @ $FFTW])
	    AC_MSG_NOTICE([Attempting linking to the user's FFTW lib ${FFTW}lib/libfftw3${float_version}])
	
	    ## if we use the user one, only link to the static lib
	    if [ test -a ${FFTW}lib/libfftw3${float_version}.a ]; then
		AC_MSG_NOTICE([Statically linking to the user's FFTW lib])
	    	CFLAGS="$CFLAGS -I${FFTW}include"
	    	LDFLAGS="$LDFLAGS ${FFTW}lib/libfftw3${float_version}.a"
		if test x$OMP_FFTW = xtrue ; then
		   if test -a ${FFTW}lib/libfftw3${float_version}_omp.a ; then
		      LDFLAGS="$LDFLAGS ${FFTW}lib/libfftw3${float_version}_omp.a"
		      AC_DEFINE([OMP_FFTW],[],[OpenMP-parallel FFT routines.])
		   fi
		fi
	        AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
	    ## Have a look for the shared object library ...
            elif [ test -a ${FFTW}lib/libfftw3${float_version}.so ]; then
		AC_MSG_NOTICE([Static link not found !! Shared object (.so) linking to the user's FFTW lib])
	    	CFLAGS="$CFLAGS -I${FFTW}include -L${FFTW}lib/"
	    	LDFLAGS="$LDFLAGS -lfftw3${float_version}"
		if test x$OMP_FFTW = xtrue ; then
		   ## check to see if it is there otherwise don't include it
		   if test -a ${FFTW}lib/libfftw3${float_version}_omp.so ; then
		      LDFLAGS="$LDFLAGS -lfftw3${float_version}_omp"
		      AC_DEFINE([OMP_FFTW],[],[OpenMP-parallel FFT routines.])
		   fi
		fi
	        AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
            else
	    	AC_MSG_ERROR([fftw3.h found, but cannot find the shared or static libraries!])
	    fi
	    ]	    
	    , [
	    ## see if we already have it	
	    AC_CHECK_HEADERS([fftw3.h],[local_fftw=true],[local_fftw=false])

	    # get desperate here and check for the header in the usual places
	    # as omp routines are special I doubt they will be installed, maybe should do a special check
	    if test x$local_fftw = "xtrue" ; then
	       LDFLAGS="$LDFLAGS -lfftw3${float_version}"
	       AC_MSG_NOTICE([Using the system's fftw])	       
	       AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
	    else
	       AC_MSG_NOTICE([No local fftw3.h found either!])
	    fi
     ]) ## closes the check headers ...
]
# closes with-FFTW
)
AM_CONDITIONAL([DEFFFTW],[$fftw])
AC_SUBST(FFTW)

######################################### GSL ###############################################
#                                                                                           #
#    GSL (GNU scientific library) can be used, link to allow for the mersenne twister RNG   #
#    and its eigenvalue stuff for large NC. Not a necessity.                                #
#                                                                                           #
#############################################################################################
gsl=false 
AC_ARG_WITH([gsl],
	AC_HELP_STRING([--with-gsl=<name>],
		       [Specify GSL location for the GF]),
            [GSL=${with_gsl}
	    ## test to see if something is there ...
	    if [ test -a ${GSL}include/gsl/gsl_rng.h ]; then
	    	AC_MSG_NOTICE([User defined Gnu Scientific Library (GSL) @ $GSL])
	  	## check if the gsl static lib is available
	    	if [ test -a ${GSL}lib/libgsl.a ]; then 
		   AC_MSG_NOTICE([Statically linking to the user's GSL lib])
	    	   CFLAGS="$CFLAGS -I${GSL}include"
	    	   LDFLAGS="$LDFLAGS ${GSL}lib/libgsl.a"	
		   AC_DEFINE([HAVE_GSL],,[Do we have the GSL libs?])
		elif [ test -a ${GSL}lib/libgsl.so ]; then
	    	   AC_MSG_NOTICE([Static link to GSL not found!! Shared object (.so) linking to the user's GSL link ])	
	    	   CFLAGS="$CFLAGS -I${GSL}include -L${GSL}lib/ "
	    	   LDFLAGS="$LDFLAGS -lgsl -lgslcblas"
		   AC_DEFINE([HAVE_GSL],,[Do we have the GSL libs?])
	    	else
	           AC_MSG_ERROR([gsl_rng.h found, but cannot find the static or shared libraries needed to link against!])
	    	fi
		echo " "
	    fi
	    ## have a look for the ones in /usr/
	    AC_CHECK_HEADERS([gsl/gsl_eigen.h],[LDFLAGS="$LDFLAGS -lgsl -lgslcblas"])
		],[ AC_MSG_NOTICE([This is a GSL-free build.]) ])
AM_CONDITIONAL([DEFGSL],[$gsl])
## Include where it is just in case people want to do something in Makefile.am
AC_SUBST(GSL)

######################################### LAPACKE ############################################
#                                                                                            #
#    Binding to the lapacke library (lapack for c), not really used too much but can be      #
#    very beneficial for large NC eigenvalues. Not a necessity.                              #
#                                                                                            #
##############################################################################################
lapacke=false 
AC_ARG_WITH([lapacke],
	AC_HELP_STRING([--with-lapacke=<name>],
		       [Specify Lapacke location for the stouting and stuff]),
            [LAPACKE=${with_lapacke}

	    if [ test -a ${LAPACKE}lapacke.h ]; then
	    	AC_MSG_NOTICE([User defined Lapacke routines])
		AC_DEFINE([HAVE_LAPACKE_H],,[Lapacke functions.])
	  	## check if the lapacke static lib is available
	    	if [ test -a ${LAPACKE}liblapacke.a ]; then 
		AC_MSG_NOTICE([Statically linking to the user's Lapacke lib])
	    	CFLAGS="$CFLAGS -I${LAPACKE}"
	    	LDFLAGS="$LDFLAGS ${LAPACKE}liblapacke.a ${LAPACKE}liblapack.a ${LAPACKE}libblas.a -lgfortran"	
	    	fi
		echo " "	   
	    else
		## have a look for the ones in /usr/
	    	AC_CHECK_HEADERS([lapacke.h],[LDFLAGS="$LDFLAGS -llapacke -llapack -lblas -lgfortran"])
	    fi
	    ],[ AC_MSG_NOTICE([This is a Lapacke-free build.]) ])
AM_CONDITIONAL([DEFLAPACKE],[$lapacke])
## Include where it is just in case people want to do something in Makefile.am
AC_SUBST(LAPACKE)

## Derivative types allowed by the code
AC_ARG_WITH([deriv],
	AC_HELP_STRING([--with-deriv=],
		       [Compile in the gauge group number NC]),
            [
	    AC_MSG_NOTICE([User specified deriv, type given as ${with_deriv}])	
	    case "${with_deriv}" in
	      deriv_full) AC_DEFINE([deriv_full],[],[Normal derivative, log def gluon fields]) ;;
	      ## JAMIE :: TODO, work out the functionals being minimised in the CG for these!
	      #deriv_fulln) AC_DEFINE([deriv_fulln],[],[Stencil derivative, log def gluon fields]) ;;
	      #deriv_linn) AC_DEFINE([deriv_linn],[],[Stencil derivative, linear gluon fields]) ;;
	      *) AC_MSG_NOTICE([Default deriv : Linear deriv with O(a) approx exponentiation (Fastest).]);;
	    esac
	    ],[ AC_MSG_NOTICE([Default deriv : Linear deriv with O(a) approx exponentiation (Fastest).]) ])

## Definitions for the fastsmearing default is the slow(stable) one
AC_ARG_ENABLE([exp_exact],
	[ --enable-exp_exact Exact exponentiation into the group for the gauge fixing],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([Exact exponentiation for the gauge transformation matrices])
       	     AC_DEFINE([exp_exact],,[Exact exponentiation for the gauge transformation matrices])   
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-exp_exact]) ;;
	esac],[])

## Definitions for the fastsmearing default is the slow(stable) one
AC_ARG_ENABLE([fastsmear],
	[ --enable-fastsmear Turn on numerically unstable fast smearing routines],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([fast and unstable routines enabled])
             AC_DEFINE([FAST_SMEAR],,[Enables faster and innaccurate smearing routines])   
	     ;;
	no) AC_MSG_NOTICE([default stable smearing routines used]) ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-fastsmear]) ;;
	esac],[])

## Finally we have routines for this node or many machines
## historically targeted at "CONDOR", default is CONDOR_MODE
AC_ARG_ENABLE([notcondor],
	[ --enable-notcondor Allows for the saving of FFT plans and some simple hashing and whatever],
	[case "${enableval}" in
	yes) ## Check the condor enable creates the dirs $prefix/Local/(Wisdom && Moments)
	     ## used for some light caching of previous results
             if [ test "$prefix" == "NONE" ] ;then
             	AC_MSG_NOTICE([Cannot save wisdom as we are not prefixed, defining CONDOR_MODE by default])
             else
		AC_MSG_NOTICE([FFTW's WISDOM routines being used, targeted specifically for this node])
		AC_DEFINE([NOT_CONDOR_MODE],,[targeted, has been seen to be quicker but do not distribute the binaries!])
		mkdir $prefix/Local/	
		mkdir $prefix/Local/Wisdom
		mkdir $prefix/Local/Moments 
             fi
	     AC_MSG_NOTICE([FFTW plans created ONCE. Subsequent transforms are read from Local/Wisdom/])
	     ;;
	no) AC_MSG_NOTICE([default condor routines called, i.e. FFTW plans created here and now]) ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-notcondor]) ;;
	esac],[])

## Include the other stuff we could possibly use
## with proj gracey
AC_ARG_WITH([GRACEY_PROJ],
	AC_HELP_STRING([--with-GRACEY_PROJ=],
		       [Compile in the specific nonexceptional projector of interest]),
            [
	    if test "x$with_GRACEY_PROJ" != "xyes" ; then	 
	       AC_MSG_NOTICE([User specified nonexceptional projector ${with_GRACEY_PROJ}.])
	       AC_DEFINE_UNQUOTED([PROJ_GRACEY], [$with_GRACEY_PROJ] , [Nonexceptional projection index.] )	
	    else
	       AC_DEFINE([PROJ_GRACEY],[0],[Nonexceptional projection index.] )	
	    fi
	    ],[ 
	    	AC_MSG_NOTICE([Default nonexceptional projector (1).])
	    	AC_DEFINE([PROJ_GRACEY],[0],[Nonexceptional projection index.] ) ] )

## Staple improvements options are IWASAKI, DBW2 , SYMANZIK and SYMANZIK_ONE_LOOP 
AC_ARG_WITH([IMPROVED_C0],
	AC_HELP_STRING([--with-IMPROVED_C0=],
		       [The 1x1 improved term])
	[
	if test x$with_IMPROVED_C0 != "xyes" ; then
	   AC_DEFINE([IWA_WEIGHT1],[$with_IMPROVED_C0],[Leading term])
	fi
	],[])

## and the rectangle terms
AC_ARG_WITH([IMPROVED_C1],
	AC_HELP_STRING([--with-IMPROVED_C1=],
		       [The 2x1 rectangle terms])
	[
	if test x$with_IMPROVED_C1 != "xyes" ; then
	   AC_DEFINE([IWA_WEIGHT1],[$with_IMPROVED_C1],[Sub-leading rectangle term])
	fi
	],[])

## if we have not specified c0 and c1
if [ test x$with_IMPROVED_C0 != "xyes" -a x$with_IMPROVED_C1 != "xyes" ] ; then
AC_ARG_WITH([IMPROVED_STAPLE],
	AC_HELP_STRING([--with-IMPROVED_STAPLE=],
		       [Compile in the type of improved staple we wish to use.]),
            [
	    AC_MSG_NOTICE([user defined ${with_IMPROVED_STAPLE} improved staple.])
	    case ${with_IMPROVED_STAPLE} in
	    DBW2) AC_DEFINE([DBW2],[],[DBW2 improved staple.]) ;;
	    IWASAKI) AC_DEFINE([IWASAKI],[],[Iwasaki improved staple.])	;;
            SYMANZIK) AC_DEFINE([SYMANZIK],[],[Symanzik improved staple.]) ;;
            *) AC_MSG_NOTICE([Unrecognised, not using rectangle staples]) ;;
	    esac
	    ],[])
fi

## are we doing overimproved smearing 
AC_ARG_WITH([OVERIMPROVED_EPSILON],
	AC_HELP_STRING([--with-OVERIMPROVED_EPSILON=],
		       [Give a value for the overimproved epsilon. Should be used in conjunction with improved staple]),
            [
	    ## if the user has left blank we do not define it, that'll teach them
	    if test "x$with_OVERIMPROVED_EPSILON" != "xyes" ; then	 
	    AC_MSG_NOTICE([User specified overimprovement factor ${with_OVERIMPROVED_EPSILON}.])
	    AC_DEFINE_UNQUOTED([epsilon],[$with_OVERIMPROVED_EPSILON], [Defined the overimproved smearing parameter epsilon.]) 
	    fi
	    ],[])

## Verbose output 
AC_ARG_ENABLE([verbiage],
	[ --enable-verbiage ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([Verbose output, mostly to look at the plaquettes at each iteration of smearing.])	
	     AC_DEFINE([verbose],[],[Verbose output printed to stdout.])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-verbiage]) ;;
	esac
	],[])

## ASCII check 
AC_ARG_ENABLE([ascii_check],
	[ --enable-ascii_check ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([If we write out an output file, this prints the output the stdout as well.])	
	     AC_DEFINE([ASCII_CHECK],[],[Ascii output to stdout.])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-ascii_check]) ;;
	esac
	],[])

## BlueGene inlined matrix multiplies
AC_ARG_ENABLE([GLU_bgq],
	[ --enable-GLU_bgq ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([Preparing for a BGQ build...])	
	     AC_DEFINE([GLU_BGQ],[],[Bluegene specific hacks.])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-GLU_bgq]) ;;
	esac
	],[])

## GIVENS APE projection
AC_ARG_ENABLE([GIVENS_APE],
	[ --enable-GIVENS_APE ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([Using the trace maximisation APE projection ...])	
	     AC_DEFINE([GIVENS_APE],[],[Trace maximisation routine for APE projection.])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-GIVENS_APE]) ;;
	esac
	],[])

## Clover improvements
AC_ARG_ENABLE([CLOVER_IMPROVE],
	[ --enable-CLOVER_IMPROVE ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([Highly improved clover term specified for topological charge measurement.])	
	     AC_DEFINE([CLOVER_IMPROVE],[],[Use the highly improved clover term.])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-CLOVER_IMPROVE]) ;;
	esac
	],[])

## and k5 to boot otherwise we default to k5 == 0.0 !!
AC_ARG_WITH([CLOVER_K5],
	AC_HELP_STRING([--with-CLOVER_K5=],
		       [Give a value for the k5 term, turns on the 3x3 staple inclusion.]),
            [ 
	    ## if the user has left blank we do not define it, that'll teach them
	    if test "x$with_CLOVER_K5" != "xyes" ; then
	    AC_MSG_NOTICE([Clover k5 = ${with_CLOVER_K5} term specified.])
	    AC_DEFINE_UNQUOTED([k5],[$with_CLOVER_K5],[K5 term in the overimproved smearing definition.]) ] 
	    fi
	    ,[])

## enable where we start our iterative topological charge measuring from
AC_ARG_WITH([TOP_VALUE],
	AC_HELP_STRING([--with-TOP_VALUE=],
		       [Tells us which measurement number from which to start checking the topological charge.]),
            [
	    ## if the user has left blank we give it a default value of 1
	    if test "x$with_TOP_VALUE" != "xyes" ; then
	    AC_MSG_NOTICE([Performing topological charge measurements from smearing iteration ${with_TOP_VALUE}.])	
	    AC_DEFINE_UNQUOTED([TOP_VALUE],[$with_TOP_VALUE],[Topological charge measurement start.]) 
	    else
	    AC_DEFINE([TOP_VALUE],[1],[Topological charge measurement start.]) 
	    fi
	    ],[])

## can now choose which RNG you want to use, boy, aren't you lucky
## options are GSL, KISS or MWC otherwise default to WELL
AC_ARG_WITH([RNG],
	AC_HELP_STRING([--with-RNG=],
		       [Specify which (p)RNG you wish to use, default is the Well.]),
            [
	    case ${with_RNG} in
	    KISS_RNG) AC_DEFINE([KISS_RNG],[],[Keep It Simple Stupid RNG (JKISS32)]) ;;
	    MWC_1038_RNG) AC_DEFINE([MWC_1038_RNG],[],[Multiply With Carry (1038) selection]) ;;
	    MWC_4096_RNG) AC_DEFINE([MWC_1038_RNG],[],[Multiply With Carry (4096) selection]) ;;
	    GSL) AC_DEFINE([GSL_RNG],[],[GSL default (Mersenne Twister) selection]) ;;
	    *) AC_MSG_NOTICE([Defaulting to the Well19937a RNG]) ;;
	    esac
	    ],[ AC_MSG_NOTICE([Defaulting to the Well19937a RNG]) ])

## enable the searching for the best gribov copies
AC_ARG_WITH([LUXURY_GAUGE],
	AC_HELP_STRING([--with-LUXURY_GAUGE=],
		       [How many random gribov copies do we look for through. Keeps the best (i.e. the one that minimises the functional).]),
            [
	    ## if the user has left blank we give it a default value of 1
	    if test "x$with_LUXURY_GAUGE" != "xyes" ; then
	    AC_MSG_NOTICE([Calculating ${LUXURY_GAUGE} Gribov copies and keeping the best.])	
	    AC_DEFINE_UNQUOTED([LUXURY_GAUGE],[$with_LUXURY_GAUGE],[Generated Gribov copies.]) 
	    else ## default to 20
	    AC_DEFINE([LUXURY_GAUGE],[20],[Generated Gribov copies.]) 
	    fi
	    ],[])

## steepest descents routines
AC_ARG_ENABLE([GLU_GFIX_SD],
	[ --enable-GLU_GFIX_SD ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([Using the SLOWER Steepest Descent gauge fixing routines])	
	     AC_DEFINE([GLU_GFIX_SD],[],[SD gauge fixing routines])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-GLU_GFIX_SD]) ;;
	esac
	],[])

## steepest descents routines
AC_ARG_ENABLE([OVERRELAX_GF],
	[ --enable-OVERRELAX_GF ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([WARNING!! Using the overrelaxed gauge fixing routines])	
	     AC_DEFINE([OVERRELAX_GF],[],[overrelaxed routines are slow])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-OVERRELAX_GF]) ;;
	esac
	],[])

## Some compiler checks
## My code is littered with consts
AC_C_CONST([])

## Also use the restrict flag a lot
AC_C_RESTRICT([])

## Many of my functions are "inlined" {compiler's choice so maybe not} ,
## hence the large binary
AC_C_INLINE([])

## Set up config.h, one of the most important headers of them all
AC_CONFIG_HEADERS([config.h])

## where are my files?
AC_CONFIG_FILES([
Makefile
src/Makefile
])

## tell us what flags we are using
echo CFLAGS used :: $CFLAGS
echo LDFLAGS used :: $LDFLAGS

AC_OUTPUT
