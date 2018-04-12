## GLU Automake instructions outline ##

## --prefix=  :: tell it where to install the binary
## --with-NC= :: number of colors, if not specified defaults to 3
## --with-ND= :: number of dimensions, if not specified defaults to 4
## --enable-single :: working precision is single, default is double
## --with-fftw= :: directory for the installed FFTW routines if --with-fftw
##                 argument is specified without a directory, will look for
##                 -lfftw3(f) in /usr/lib/ 
## --with-gsl= :: directory for GSL
## --with-lapacke= :: directory for lapacke install, used for large-NC eigen
##                    operations, good binding with complex.h
## --enable-notcondor :: enable the single-machine cache-ing and architecture
##                       dependent optimisations
## --with-GRACEY_PROJ= :: MOMggg scheme projection number, 0 is default
## --with-TOP_VALUE= :: Specify number of smearing steps before measuring the
##                      topological charge 
## --enable-GIVENS_APE :: use the Cabbibo-Marinari trace-maximisation for
##                        the APE and HYP smearings
## --enable-GLU_bgq :: inlined matrix multiplies via macros
## --with-IMPROVED_C0=  :: leading 1x1 term if not using IWASAKI, SYMANZIK or DBW2
## --with-IMPROVED_C1=  :: subleading 2x1 term only effective if C0 is defined
## --with-IMPROVED_STAPLE= :: Rectangle coefficients, choices are :
## 			     IWASAKI, SYMANZIK and DBW2
## --with-OVERIMPROVED_EPSILON= :: specify a value for the overimprovement
##                                 factor, default is 0.0
## --enable-verbiage :: more verbose output than usual
## --enable-ascii_check :: prints to stdout the data we would write to our
##                         binary files
## --enable-PLAQUETTE_FMUNU :: plaquette defintition of field strength tensor
## --enable-CLOVER_IMPROVE :: Are we using the O(a^4) clover def?
## --with-CLOVER_K5= :: specify the multiplicative term for the 3x3 wilson
##                      loop for the heavily improved field strength tensor
## --with-RNG= :: choices are KISS, MWC and GSL. Default is the WELL
## --with-LUXURY_GAUGE= :: choose the "best" (one that minimises the
## 		       	  functional) from specified number of random
##                         initial transformations, default is 20
## --enable-GLU_GFIX_SD :: turn on the steepest descent codes (default is the CG)
## --enable-OVERRELAXED_GF :: switch to the overrelaxed codes for gauge fixing

AC_INIT([GLU], [1.1], [renwick.james.hudspith@gmail.com])
AM_INIT_AUTOMAKE([-Wall subdir-objects])
AC_PROG_CC
AM_PROG_CC_C_O
AC_LANG([C]) 
AC_PROG_RANLIB ## include a check for libtool if not use ranlib?
AM_PROG_AR
AC_C_BIGENDIAN 

## SET up the m4
AC_CONFIG_MACRO_DIR([m4])

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
    AM_CONDITIONAL([DOX],[false])
    AM_CONDITIONAL([DOT],[false])
    AM_CONDITIONAL([PDFLTX],[false])
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

## Compile for a set NC, default is NC=3
AC_ARG_WITH([NC],
	AC_HELP_STRING([--with-NC=],
		       [Compile in the gauge group number NC]),
            [
	    AC_MSG_NOTICE([User specified NC, compiling for SU(${with_NC})])
	    AC_DEFINE_UNQUOTED([NC], [${with_NC}] , [Compiled for SU(NC)] )	
	    ],[
	    AC_MSG_NOTICE([User unspecified NC, default to SU(3)])	
	    ])

## Compile for a set ND, default is ND=4
AC_ARG_WITH([ND],
	AC_HELP_STRING([--with-ND=],
		       [Compile in the number of dimensions ND]),
            [
	    AC_MSG_NOTICE([User specified ND, compiling for ND=${with_ND}])
	    AC_DEFINE_UNQUOTED([ND], [${with_ND}] , [Compiled for ND] )	
	    ],[
	    AC_MSG_NOTICE([User unspecified ND, default to 4])	
	    ])

## C standard libraries
AC_CHECK_HEADERS([math.h],[],AC_MSG_ERROR([We implicitly use math.h functions and cannot find it!]))

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

## check for immintrin.h -> Intel vector instructions
AC_CHECK_HEADERS([immintrin.h],
	[AC_MSG_NOTICE([Vector intrinsics being used])],[])

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
##                                                                                           ##
##    This one is tricky, looks first for the static link to FFTW library, if not found      ##
##    then the shared. If --with-fftw= then it looks in /usr/lib for the shared              ##
##    otherwise we miss out on all of the Fourier transform stuff like Fourier acceleration, ##
##    Momentum space gluons and quenched U(1). Not a necessity.                              ##
##                                                                                           ##
###############################################################################################

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

	    ## get desperate here and check for the header in the usual places
	    ## as omp routines are special I doubt they will be installed, maybe should do a special check
	    if test x$local_fftw = "xtrue" ; then
	       LDFLAGS="$LDFLAGS -lfftw3${float_version}"
	       AC_MSG_NOTICE([Using the system's fftw])	       
	       AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
	    else
	       AC_MSG_NOTICE([No local fftw3.h found either!])
	    fi
     ]) ## closes the check headers ...
]
## closes with-FFTW
)
AM_CONDITIONAL([DEFFFTW],[$fftw])
AC_SUBST(FFTW)

######################################### GSL ################################################
##                                                                                          ##
##    GSL (GNU scientific library) can be used, link to allow for the mersenne twister RNG  ##
##    and its eigenvalue stuff for large NC. Not a necessity.                               ##
##                                                                                          ##
##############################################################################################
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
##                                                                                          ##
##    Binding to the lapacke library (lapack for c), not really used too much but can be    ##
##    very beneficial for large NC eigenvalues. Not a necessity.                            ##
##                                                                                          ##
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

## Use exact exponential in the GF
AC_ARG_ENABLE([exp_exact],
	[ --enable-exp_exact Exact exponentiation into the group for the gauge fixing],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([Exact exponentiation for the gauge transformation matrices])
       	     AC_DEFINE([exp_exact],,[Exact exponentiation for the gauge transformation matrices])   
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-exp_exact]) ;;
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
	    if [ test "x$with_GRACEY_PROJ" != "xyes" ] ; then	 
	       AC_MSG_NOTICE([User specified nonexceptional projector ${with_GRACEY_PROJ}.])
	       AC_DEFINE_UNQUOTED([PROJ_GRACEY], [${with_GRACEY_PROJ}] , [Nonexceptional projection index.] )	
	    else
	       AC_DEFINE([PROJ_GRACEY],[0],[Nonexceptional projection index.] )	
	    fi
	    ],[ 
	    	AC_MSG_NOTICE([Default nonexceptional projector (1).])
	    	AC_DEFINE([PROJ_GRACEY],[0],[Nonexceptional projection index.] ) ] )

## c0 is the term in front of the elementary staple
AC_ARG_WITH([IMPROVED_C0],
	AC_HELP_STRING([--with-IMPROVED_C0=],
		       [The 1x1 improved term]),
	[
	if [ test "x$with_IMPROVED_C0" != "xyes" ] ; then
	   AC_DEFINE_UNQUOTED([IWA_WEIGHT1],[${with_IMPROVED_C0}],[Leading term])
	fi
	## otherwise do nothing
	],[])

## and the rectangle terms
AC_ARG_WITH([IMPROVED_C1],
	AC_HELP_STRING([--with-IMPROVED_C1=],
		       [The 2x1 rectangle terms]),
	[
	if [ test "x$with_IMPROVED_C1" != "xyes" ] ; then
	   AC_DEFINE_UNQUOTED([IWA_WEIGHT2],[${with_IMPROVED_C1}],[Sub-leading rectangle term])
	fi
	],[])

## if we have not specified c0 and c1
if [ test "x$with_IMPROVED_C0" != "xyes" -a "x$with_IMPROVED_C1" != "xyes" ] ; then
AC_ARG_WITH([IMPROVED_STAPLE],
	AC_HELP_STRING([--with-IMPROVED_STAPLE=],
		       [Compile in the type of improved staple we wish to use.]),
            [
	    AC_MSG_NOTICE([user defined "${with_IMPROVED_STAPLE}" improved staple.])
	    case "${with_IMPROVED_STAPLE}" in
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
	    if [ test "x$with_OVERIMPROVED_EPSILON" != "xyes" ] ; then	 
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

## plaquette field strength tensor
AC_ARG_ENABLE([PLAQUETTE_FMUNU],
	[ --enable-PLAQUETTE_FMUNU ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([Simple plaquette definition of the field strength tensor.])	
	     AC_DEFINE([PLAQUETTE_FMUNU],[],[Plaquette field strength tensor.])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-PLAQUETTE_FMUNU]) ;;
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
	    MWC_4096_RNG) AC_DEFINE([MWC_4096_RNG],[],[Multiply With Carry (4096) selection]) ;;
	    XOR_1024_RNG) AC_DEFINE([XOR_1024_RNG],[],[XOR (1024) selection]) ;;
	    GSL) AC_DEFINE([GSL_RNG],[],[GSL default (Mersenne Twister) selection]) ;;
	    WELL_512_RNG) AC_DEFINE([WELL_512_RNG],[],[Well (512) rng selection]) ;;
	    *) AC_MSG_NOTICE([Defaulting to the MWC_4096_RNG]) ;;
	    esac
	    ],[ AC_MSG_NOTICE([Defaulting to the MWC_4096_RNG]) ])

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

## enable the searching for the best gribov copies
AC_ARG_WITH([NBLOCK],
	AC_HELP_STRING([--with-NBLOCK=],
		       [Innermost loop unrolling factor for matrix multiplies.]),
            [
	    ## if the user has left blank we give it a default value of 1
	    if test "$with_NC" -gt 3 ; then
	       if test "x$with_NBLOCK" != "xyes" ; then
	       	  if test $with_NBLOCK -gt $with_NC ; then
		     AC_MSG_ERROR([NBLOCK of ${with_NBLOCK} cannot be greater than SU($with_NC).])
		  else
		     let 'check_odd=with_NC%2'
	       	     if test $check_odd = 1 ; then
	                let 'check=(with_NC-1)%with_NBLOCK'
	       	     else
		        let 'check=with_NC%with_NBLOCK'
	       	     fi  
	       	     if test $check != 0 ; then
	                AC_MSG_ERROR([NBLOCK of ${with_NBLOCK} not a factor of SU($with_NC).])
	       	     else
		        AC_MSG_NOTICE([Performing ${with_NBLOCK} inner loop unrolling.])
		        AC_DEFINE_UNQUOTED([NBLOCK],[$with_NBLOCK],[Loop unrolling factor.])	
	       	    fi
		 fi
	       else ## default to 1 which is no loop unrolling
	       	  AC_DEFINE([NBLOCK],[1],[No loop unrolling.]) 
	       fi
	    else
		if test "x$with_NBLOCK" != "xyes" ; then
	       	   AC_MSG_NOTICE([Using default inline matrix routines. Ignoring blocking factor.])
		fi
	    fi	    
	    ],[])

## only timing the gradient flow
AC_ARG_ENABLE([WFLOW_TIME_ONLY],
	[ --enable-WFLOW_TIME_ONLY ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([Only timing the gradient flow])	
	     AC_DEFINE([WFLOW_TIME_ONLY],[],[Only timing the gradient flow])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-WFLOW_TIME_ONLY]) ;;
	esac
	],[])
	

## what time to stop the wflow
AC_ARG_WITH([WFLOW_TIME_STOP],
	AC_HELP_STRING([--with-WFLOW_TIME_STOP=],
		       [Specify at what time you want the flow to stop]),
            [
	    if test "x$with_WFLOW_TIME_STOP" != "xyes" ; then
	        AC_MSG_NOTICE([Stopping flow at t=${with_WFLOW_TIME_STOP}.])
	        AC_DEFINE_UNQUOTED([WFLOW_TIME_STOP],[$with_WFLOW_TIME_STOP],[Wilson flow maximum time stopping condition.])
	    else
		AC_MSG_NOTICE([If you want to use WFLOW_TIME_STOP you should specify a stopping time.])
	    fi
	    ],[])

## what time to stop the wflow
AC_ARG_WITH([WFLOW_MEAS_START],
	AC_HELP_STRING([--with-WFLOW_MEAS_START=],
		       [Specify at what time you want to start measuring flow observables]),
            [
	    if test "x$with_WFLOW_MEAS_START" != "xyes" ; then
	        AC_MSG_NOTICE([Measuring from t=${with_WFLOW_MEAS_START}.])
	        AC_DEFINE_UNQUOTED([WFLOW_MEAS_START],[$with_WFLOW_MEAS_START],[Wflow observable measurements start from here])
	    else
		AC_MSG_NOTICE([If you want to use WFLOW_MEAS_START you should specify a starting time.])
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

## steepest descents OverRelaxed routines
AC_ARG_ENABLE([OVERRELAXED_GF],
	[ --enable-OVERRELAXED_GF ],
	[case "${enableval}" in
	yes) AC_MSG_NOTICE([WARNING!! Using the overrelaxed gauge fixing routines])	
	     AC_DEFINE([OVERRELAXED_GF],[],[overrelaxed routines are slow])
	     ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-OVERRELAXED_GF]) ;;
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
tests/Makefile
])

## tell us what flags we are using
echo CFLAGS used :: $CFLAGS
echo LDFLAGS used :: $LDFLAGS

AC_OUTPUT
