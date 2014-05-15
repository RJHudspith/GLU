#ifndef CONFIG_H
#define CONFIG_H

//#pragma clang diagnostic ignored "-Wunknown-pragmas"

// blank here as this is now being pushed over to automake
#define HAVE_PREFIX "/home/jamie/PHYSICS/GaugeFixer/"
#define NOT_CONDOR_MODE
//#define HAVE_CILK_H superseded by openmp pragmas
//#define OMP_FFTW
#define HAVE_OMP_H

//#define HAVE_GSL
//#define HAVE_LAPACKE_H
#define NC 18
#define ND 4
//#define ASCII_CHECK

// define for fftw.h
#define HAVE_FFTW3_H

//#define GLU_BGQ
//#define verbose

//#define SINGLE_PREC

#define HAVE_UNISTD_H
#define HAVE_STDINT_H
#define HAVE_SYS_TIME_H
#define HAVE_TIME_H

//#define GLU_GFIX_SD

//#define U1_DHT

//#define CHROMA_APE just included as a check
//#define GIVENS_APE

// can define WRITE_FIELDS or EXCEPTIONAL, else will compute the non-exceptional
//g2 and g3 procedure...

#define PSQ_MOM
//#define SIN_MOM
#define PROJ_GRACEY 0

#define WORDS_BIGENDIAN 0

//#define FAST_SMEAR sometimes this is useful

// choices are CLOVER_IMPROVE or not ...
// should put these in gftests and clover ...
// -- should be options in the input file, slow routines should not be guarded
//#define CLOVER_IMPROVE
//#define PLAQUETTE_FMUNU
//#define k5 0.0

//#define TOP_VALUE 1

// both of these should be taken from the input file? Maybe not, only eta?
//#define SYMANZIK
//definition of the eta parameter
//#define epsilon -1.0

//#define MWC_4096_RNG
//#define MWC_1038_RNG
//#define KISS_RNG

// -- handled as enables in the automake
//gtrans defines
//#define LUXURY_GAUGE 40
//#define CAREFUL 1

//#define deriv_linn
//#define deriv_fullnn
//#define deriv_full

//#define nn1 4.0/3.0
//#define nn2 -1.0/12.0

//#define exp_a2_approx
//#define exp_exact // cannot use in conjunction with MLG !!!

//#define IWA_WEIGHT1 1.0
//#define IWA_WEIGHT2 -1.0

#endif
