/*
    Copyright 2013 Renwick James Hudspith

    This file (GLU_definitions.h) is part of GLU.

    GLU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GLU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GLU.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   @file GLU_definitions.h 
   @brief I needed somewhere to store all of the definitions I was using
 **/

#ifndef GLU_DEFINITIONS_H
#define GLU_DEFINITIONS_H

#include "Mainfile.h"

/**
   @param OMP_FFTW
   @brief uses the openmp'd FFT routines
   Instead of performing the element by element naive parallelism that
   I originally implemented this uses FFTW's own openMP'd routines
   This is of great benefit at large volumes!
 */
#ifdef OMP_FFTW
  #ifndef HAVE_OMP_H
    #undef OMP_FFTW
  #endif
#endif

#ifdef SINGLE_PREC
/**
   @typedef GLU_complex
   @brief single precision complex variables
   @typedef GLU_real
   @brief single precision variables
 */
   typedef float complex GLU_complex ;
   typedef float GLU_real ;
   #define fftw_free fftwf_free
   #define fftw_destroy_plan fftwf_destroy_plan   
   #define fftw_cleanup fftwf_cleanup
   #define fftw_malloc fftwf_malloc
   #define fftw_plan fftwf_plan
   #define fftw_execute fftwf_execute
   #define fftw_init_threads fftwf_init_threads
   #define fftw_cleanup_threads fftwf_cleanup_threads
#ifdef OMP_FFTW
   #define fftw_plan_with_nthreads fftwf_plan_with_nthreads 
#endif
   #define fftw_import_wisdom_from_file fftwf_import_wisdom_from_file
   #define fftw_export_wisdom_to_file fftwf_export_wisdom_to_file
   #define fftw_plan_dft fftwf_plan_dft
   #define fftw_plan_r2r fftwf_plan_r2r
   #define fftw_r2r_kind fftwf_r2r_kind
   #define creal crealf
   #define cimag cimagf
   #define cabs cabsf
#else
/**
   @typedef GLU_complex
   @brief double precision complex variables
   @typedef GLU_real
   @brief double precision variables
 */
   typedef double complex GLU_complex ;
   typedef double GLU_real ;
#endif

/**
   @def GLU_RESTRICT
   @brief restrict flag
 */
#define GLU_RESTRICT __restrict

/**
   defines for our control
   @def GLU_FAILURE
   @brief flag for if something did not work
   @def GLU_SUCCESS
   @brief flag for if something worked as expected
 */
#define GLU_FAILURE -1
#define GLU_SUCCESS !GLU_FAILURE

// some generic constants I use all over the place
#define OneOI2 -0.5 * I
#define OneO3 0.3333333333333333

/**
   @def MPI
   @brief Pi from math.h
   @def TWOPI 
   @brief is \f$ 2.0\times\pi \f$
   @def PIOtwo
   @brief is \f$ 0.5\times\pi \f$
 */
#define MPI M_PI //3.141592653589793
#define TWOPI 6.283185307179586
#define PIOtwo M_PI_2 //1.5707963267948966

/**
   @def r2
   @brief the \f$ \sqrt{2} \f$ in GLU_real precision
 */
#define r2 M_SQRT2 // sqrt( 2. )

/**
   sin tolerances for the exact exponentiation ..
   @def SINTOL 
   @brief sin tolerance for the exact exponentiation
   @def STOL
   @brief sin tolerance for the log exact exponentiation
   @def SINTOLSU2
   @brief we can be more lax with SU(2), in fact we usually have to be
 */
#define SINTOL 0.05
#ifdef SINGLE_PREC
  #define STOL 1E-6
  #define SINTOLSU2 0.05 // this is for the su(2)
#else
  #define STOL 1E-13
  #define SINTOLSU2 1E-4 // this is for the su(2)
#endif

/**
   @def DBL_MIN
   @brief the smallest representable double used forexact exponentiation tolerance
 */
#ifndef __DBL_MIN__
  #define DBL_MIN 2.2250738585072014e-308
#else
  #define DBL_MIN __DBL_MIN__
#endif

/**
   @def FLT_MIN
   @brief the smallest representable float
 */
#ifndef __FLT_MIN__
  #define FLT_MIN 1.17549435e-38F
#else
  #define FLT_MIN __FLT_MIN__
#endif

/**********************************

       Lattice geometry types

***********************************/

/**
   @def NC
   @brief number of colors of our theory, if not specified default to 3
 */
#ifndef NC
  #define NC 3
#endif

/**
   @def NCNC
   @brief we only use square matrices, so they are always \f$ NC^2 \f$
   size of our gauge matrices
 */
#define NCNC NC * NC

/**
   @def ND
   @brief dimensions of our theory, if not specified default to 4
 */ 
#ifndef ND
  #define ND 4
#endif

/**
   @def HERMSIZE
   @brief upper triangular of a hermitian matrix plus the NC-1 diagonal terms 
   is enough to describe all the data

   this uses rounding down of integer division.
 */
#define HERMSIZE (int)( ( NC * ( NC + 1 ) >> 1 ) - 1 )

/**
   @def TRUE_HERM
   @brief cheat for SU(3) speedup
   @warning brackets are a must!
 */
#if NC == 3
#define TRUE_HERM ( HERMSIZE-1 )
#else
  #define TRUE_HERM HERMSIZE
#endif

/**
   @def INLINE_VOID
   @brief tell the compiler to inline specific functions if they have been loop unrolled by hand 
   @def INLINE_STATIC_VOID
   @brief synactical sugar 
   @def INLINE_DOUBLE_COMPLEX
   @brief synactical sugar
   only inline stuff I have loop unrolled... NC < 4 at the moment
 */
#if NC < 4
  #define INLINE_VOID inline void
  #define INLINE_STATIC_VOID inline static void
  #define INLINE_GLU_COMPLEX inline GLU_complex 
#else
  #define INLINE_VOID void
  #define INLINE_STATIC_VOID static void
  #define INLINE_GLU_COMPLEX GLU_complex 
#endif

/**
   @def HAVE_CILK_H
   @brief if we want to use the cilk shared memory parallelism we create these
   names.
   @warning I do not have the same control with reductions as with openMP
   @def PFOR
   @brief CILK parallel form, or not
   @def PSPAWN
   @brief CILK spawns a thread, or not
   @def PSYNC
   @brief CILK collects the threads together, or not

   set up a parallel for, and spawn if we have the cilk library
   This is essentially redundant as we now have OpenMP bindings.
   need to control this behaviour if using OMP ...
 */
#ifdef HAVE_CILK_H
  #define PFOR cilk_for
  #define PSPAWN cilk_spawn
  #define PSYNC cilk_sync
#else
  #define PFOR for
  #define PSPAWN
  #define PSYNC
#endif

/**
   @def LVOLUME
   @brief macro-ises the lattice volume
   set up the lattice volume ...
 */
#ifdef VOLUME
  #define LVOLUME VOLUME
#else
  #define LSQ Latt.Lsq
  #define LCU Latt.Lcu
  #define LVOLUME Latt.Volume
#endif 

/******************************************************

  Stuff for the header of the NERSC config ( IO/{}.c )

*******************************************************/

/**
   @def STRING_SIZE
   @brief sets the string size for the input reader (input_reader.c)
 */
#define STRING_SIZE 264 

/**
   @def PLAQ_AND_TRACE_TOL
   @brief tolerance before complaint for reading in the NERSC header
   I output this to a much higher precision but if the
   configs being read are floats we don't have much choice
 */
#define PLAQ_AND_TRACE_TOL 1E-6

/**
   @def MAX_LINE_LENGTH
   @brief 128 bytes is quite long for a file
   @def MAX_TOKENS
   @brief maximum number of tokens in the header I only output 21
   stuff from chklat.c for reading the QCDheader 
 */
#ifndef MAX_LINE_LENGTH
  #define MAX_LINE_LENGTH 128
#endif
#define MAX_TOKENS 36

/**
   @var typedef unsigned int uint32_t
   @brief define the type uint32_t

   chklat.c defines if for some odd reason there is no stdint.h
 */
#ifndef HAVE_STDINT_H
  #ifndef __INT_MAX__
    #define INT_MAX 2147483647 // This is its value in limits.h maximum signed int
  #else
    #define INT_MAX __INT_MAX__
  #endif
  #if (INT_MAX==2147483647)
    #define INTS_ARE_32BIT
    typedef unsigned int uint32_t ;
  #else
    #undef INTS_ARE_32BIT
    typedef unsigned short uint32_t ;
  #endif
typedef unsigned long int uint64_t ;
#endif

/*************************************

  Initial set-up and definitions

**************************************/

/**
   @def NOT_CONDOR_MODE
   @brief if we are not running embarassingly parallel

   Are we in "NOT_CONDOR_MODE" where we use the FFTW wisdom and other 
   information-caching tactics to speed up the procedure.

   If we are not distributing jobs to many different machines we can compile
   in some timesavers for this machine. These include saving FFTW plans.
   Saving momentum lists and suchlike.
 */
#ifndef NOT_CONDOR_MODE
  #define CONDOR_MODE
#endif

/**
   @def WORDS_BIGENDIAN
   @brief syntactical sugar for endian-ness
   Configure script sets WORDS_BIGENDIAN to 1 if big_endian otherwise
   doesn't define it, we turn it into a more true/false thing here
 */
#ifndef WORDS_BIGENDIAN
  #define WORDS_BIGENDIAN 0
#endif

/**
   @def NSU2SUBGROUPS
   @brief number of su(2) subgroup embeddings
 */
#define NSU2SUBGROUPS (NC*(NC-1))/2

/****************************************

   Cut routines and the like (Cuts/{}.c )

*****************************************/

/**
   @def SIN_MOM
   @brief use the 2sin(p/2) def of the momentum for projectors and whatnot
*/
#ifndef PSQ_MOM
  #define SIN_MOM
#endif

/**
   @def PREC_TOL
   @brief precision tolerance from zero of the Gracey projector ...
*/
#ifdef SINGLE_PREC
  #define PREC_TOL NC*(1.0E-6)
#else
  #define PREC_TOL NC*(1.0E-14)
#endif

/**********************************

     Defines for Landau/gtrans.c 

 **********************************/
 
/**
   @def deriv_fullnn
   @brief log_def of the fields for gauge fixing with a next-neares stencil improved derivative
   @def deriv_fulln
   @brief log_def of the fields for gauge fixing with a stencil derivative
   @def deriv_full
   @brief log_def of the fields for gauge fixing with a normal derivative
   @def deriv_linn
   @brief log_def of the fields for gauge fixing with a stencil derivative
   @def deriv_lin
   @brief log_def of the fields for gauge fixing with a normal derivative
   @def deriv_MLG
   @brief for SU(2) we can stereographically project our fields to create the
   maximal landau gauge (MLG).

   default to the usual fast AntiHermitian_proj deriv/
 */
#if !( defined deriv_fulln) && !(defined deriv_full) && !(defined deriv_linn) && !(defined deriv_fullnn) 
  // lin def is the default behaviour
  #define deriv_lin 
#endif

#if ( defined deriv_full ) || ( defined deriv_fulln )
  // default behaviour ... hmm, if we really want we should define APPROX_LOG 
  // somewhere else
  #define WARMUP_DERIV
#endif

/**
   defs for the exponentiation used in the steepest descents
   @def exp_exact 
   @brief exact exponentiation in the gauge fixing
   @def exp_a2_approx
   @brief adding the extra A^2 exponential term and reunitarising
   @def exp_approx
   @brief expanding the exponential to the first term an reunitarising
 */
#if !( defined exp_exact ) && !( defined  exp_a2_approx )
    #define exp_approx // the usual 1 + da term with reunitarisation
#endif

/******************************************************

Defines for the gauge fixing routines ( Landau/{}.c )

******************************************************/

/**
   @def SLOW_GF
   @brief default is that this is off. Uses slightly less memory only used in Landau.h
 */
#ifndef SLOW_GF
  #define FAST_GF
#endif

/**
   @def FA
   @brief default is that this is off and FA is defined. Steepest descent gauge fixing used in Coulomb.h and Landau.h
 */
#ifndef SD
  #define FA
#endif

/**
   @def SMACC
   @brief slightly smaller GF accuracy for the smearing-prec because it doesn't matter that much used in Landau.h
 */
#define SMACC 1E-12

/**
   @def GF_GLU_FAILURES
   @brief number of re-attempts to fix the gauge before complaint. used in functions in Coulomb.h and Landau.h
 */
#ifndef GF_GLU_FAILURES
  #define GF_GLU_FAILURES 7
#endif

/**
   @def MAX_LANDAU
   @brief \f$ p^2 \f$ max in our FA code
   @def MAX_COULOMB
   @brief \f$ p^2 \f$ max for one dim less

   The factor of ND in #MAX_LANDAU and #MAX_COULOMB comes from from our measurement of \f$ p^2 \f$ <br>
   \f[

       2.0 ( ND - \sum_{i=0}^{ND-1} \cos\left( \frac{p_i \: \pi}{L_i}\right) )

       \f]
 */
#define MAX_LANDAU ( 4. * ND / (double)LVOLUME )
#define MAX_COULOMB ( 4. * ( ND - 1 ) / (double)LCU ) 

/**
   @def GFNORM_LANDAU
   @brief normalisation factor for the GF accuracy
   @def GFNORM_COULOMB
   @brief normalisation factor for the GF accuracy

   normalisation factors for our gauge fixing, factor of two from the cheating used in the calculation of the trace same in #GFNORM_LANDAU #GF_COULOMB
 */
#define GFNORM_LANDAU ( 1. / (double)( NC * LVOLUME ) )
#define GFNORM_COULOMB ( 1. / (double)( NC * LCU ) )

/**
   @def nn1 
   @breief correction term for the standard finite difference
   @def nn2
   @breief correction term for the finite difference from the next nearest neighbour
 */
#ifndef nn1
  #define nn1 1.125     // BOWMAN == 4.0 / 3.0 
#endif
#ifndef nn2
  #define nn2 -1.0/24.0 // BOWMAN == -1.0/12.0
#endif

/**
   @def nnn1 
   @breief correction term for the standard finite difference
   @def nnn2
   @breief correction term for the finite difference from the next nearest neighbour
   @def nnn3
   @breief correction term for the finite difference from the next-next nearest neighbour
 */
#define nnn1 75.0/64.0 // BOWMAN == 49./36.0
#define nnn2 -25.0/384.0 // BOWMAN == -5.0/36.0
#define nnn3 3.0/640.0 // BOWMAN == 1.0/90.0

/**
   @def WORST_COPY
   @brief the worst gribov copy, from its gauge functional
   @def BEST_COPY
   @brief the best gribov copy, from its gauge functional
 */
#ifdef LUXURY_GAUGE
 #ifndef WORST_COPY
  #define BEST_COPY
 #endif
#endif

/**
   @def CG_MAXITERS
   @brief Conjugate Gradient Maximum Iterations before an Steepest Descent step
 */
#ifndef CG_MAXITERS
 #define CG_MAXITERS 35
#endif

/**
   @def CG_TOL
   @brief specific tolerances for the CG
 */
#ifdef SINGLE_PREC
  #define CG_TOL 1E-6
#else
  #define CG_TOL 5E-12
#endif

/***********************************************

  Defines for the smearing routines (Smear/{}.c)

************************************************/

// unless we explicitly want the givens rotation APE smearing
#ifndef GIVENS_APE
  #define N_APE //rescaled projection is the best projection 
#endif

// define this because we are using iterative methods in the exponential
#if !( defined HAVE_LAPACKE_H || defined HAVE_GSL_H )
  #define SCHED schedule(dynamic)
#else // leave it blank
  #define SCHED 
#endif

/**
   @def one_min_a1
   @brief either ( 1 - Latt.sm_alpha[0] ) which is the standard, or just 1.0. Depends on whether CHROMA_APE is defined
   @def one_min_a2
   @brief either ( 1 - Latt.sm_alpha[1] ) which is the standard, or just 1.0. Depends on whether CHROMA_APE is defined
   @def one_min_a3
   @brief either ( 1 - Latt.sm_alpha[2] ) which is the standard, or just 1.0. Depends on whether CHROMA_APE is defined

   @def alpha1
   @brief the "normalised" smearing parameter for the links
   @def alpha2
   @brief the "normalised" smearing parameter for the level1 links
   @def alpha3
   @brief the "normalised" smearing parameter for the level2 links

   Chroma uses a weird projection ...
   convert our smearing alpha's to their "perturbative form"
 */
#if ND == 4
// chroma's conditions for APE are weird/closer to the original APE paper
 #ifdef CHROMA_APE
  #define one_min_a1 Latt.sm_alpha[0]
  #define one_min_a2 Latt.sm_alpha[1]
  #define one_min_a3 Latt.sm_alpha[2]
  #define alpha1 1.0
  #define alpha2 1.0
  #define alpha3 1.0
#else 
  #define one_min_a1 ( 1.0 - Latt.sm_alpha[0] )
  #define one_min_a2 ( 1.0 - Latt.sm_alpha[1] )
  #define one_min_a3 ( 1.0 - Latt.sm_alpha[2] )
  #define alpha1 Latt.sm_alpha[0] / 6.0
  #define alpha2 Latt.sm_alpha[1] * 0.25
  #define alpha3 Latt.sm_alpha[2] * 0.5
 #endif
#endif

/**
   @def SLOW_SMEAR
   @brief the default, does not use the approximate exponentiation or anything untoward.
 */
#ifndef FAST_SMEAR
  #define SLOW_SMEAR
#endif

/**
   @def epsilon
   @brief eta is the overimprovement factor in the improved smearing routines
   found in <a href="linkURL">http://arxiv.org/abs/0801.1165 </a>
   if it is not defined set to usual improvement factor of 0.0
 */
#ifndef epsilon
  #define epsilon 0.0
#endif

/**
   @def IMPROVED_SMEARING
   @brief if we are using the (over) improved smearing routines this is defined
   @def IWA_WEIGHT1
   @brief the first weight of the "improved smearing" multiplies the tree staples
   @def IWA_WEIGHT2
   @brief the second weight of the "improved smearing" multiplies the corrective rectangles

   over-improvement types available are IWASAKI, DBW2, SYMANZIK, SYMANZIK_ONE_LOOP 
 */
#if ( defined IWA_WEIGHT1 ) && ( defined IWA_WEIGHT2 )
  #define IMPROVED_SMEARING
#else
  #ifdef SYMANZIK
    #define IMPROVED_SMEARING
    #define IWA_WEIGHT1 ( 5.0 - 2.0 * epsilon ) / 3.
    #define IWA_WEIGHT2 -( 1.0 - epsilon ) / 12.  
  #elif defined IWASAKI
    #define IMPROVED_SMEARING
    #define IWA_WEIGHT1 ( 1.0 + 2.648 * ( 1.0 - epsilon ) )
    #define IWA_WEIGHT2 ( -0.331 * ( 1.0 - epsilon ) )
  #elif defined DBW2
    #define IMPROVED_SMEARING
    #define IWA_WEIGHT1 ( 1.0 + 11.2536 * ( 1.0 - epsilon ) )
    #define IWA_WEIGHT2 ( -1.4069 * ( 1.0 - epsilon ) )
  #elif defined SYMANZIK_ONE_LOOP
    #define IMPROVED_SMEARING
    extern GLU_real improve ;
    #define ONELOOPSYM_FACT -0.32590381274870533
    #define ALPHAS ONELOOPSYM_FACT * log( improve )
    #define SYMONELOOP_WEIGHT1 1.
    #define SYMONELOOP_WEIGHT2 -( 1.0 - epsilon ) * ( 1. + 0.485 * ALPHAS ) / ( 20. * sqrt( improve ) ) 
    #define SYMONELOOP_WEIGHT3 -( 1.0 - epsilon ) * 0.03325 * ALPHAS / sqrt( improve )
#endif
#endif

/*******************************************************
  O(a^4) clover-improvement measure ( Landau/clover.c )
********************************************************/

/**
   @def NK5
   @brief if k5 is 0.0 which is the standard, we use the shorter clover factors
 */
#ifndef k5
  #define NK5
#endif

/**
   @def Clover_k5
   @brief the term that multiplies the 3x3 wilson loop
   @def Clover_k1
   @brief the term that multiplies the 1x1 wilson loop
   @def Clover_k2
   @brief the term that multiplies the average of the 1x2 wilson loops
   @def Clover_k3
   @brief the term that multiplies the average of the 1x3 wilson loops
   @def Clover_k4
   @brief the term that multiplies the 2x2 wilson loop

   All of the factors #Clover_k5 #Clover_k4 #Clover_k3 #Clover_k2 #Clover_k1 are used in the improved clover definition in clover.c
 */
// factors for the overimprovement, global variables used if tadpole improving
#ifdef NK5
  #define Clover_k5 0. // just keep this at 0.
  #define Clover_k1 ( 19./9. ) 
  #define Clover_k2 ( 1./36. )
  #define Clover_k3 0.5 * ( -32./45. ) // half for the 1x2,2x1 average
  #define Clover_k4 0.5 * ( 1./15. ) // half for the 1x3,3x1 average
#else
  #define Clover_k1 ( 19./9. - ( 55. * k5 ) )
  #define Clover_k2 ( 1./36. - ( 16. * k5 ) )
  #define Clover_k3 ( ( ( 64. * k5 ) - 32./45. ) * 0.5 )
  #define Clover_k4 ( ( 1./15. - ( 6. * k5 ) ) * 0.5 )
#endif

/**
   @def CL_TOL
   @brief if the k{n} term for the highly improved clover from clover.h is less than this we ignore it
*/
#define CL_TOL 1E-8

/***************************

  Control for the output

***************************/

/**
   @def OUT_BIG
   @brief define the endianess to write out to
   can either be OUT_BIG or OUT_LITTLE
 */
#ifndef OUT_LITTLE
  #define OUT_BIG
#endif

/**
   @def OUT_DOUBLE
   @brief define the precision we write our configs in
   must be either OUT_SINGLE or OUT_DOUBLE
 */
#ifndef OUT_SINGLE
 #define OUT_DOUBLE
#endif

/**************************

    Control for the RNG

**************************/

/**
   @def UINT_MAX
   @brief maximum unsigned integer
 */
#ifndef __UINT_MAX__
  #define UINT_MAX 0xffffffff // maximum unsigned int
#else
  #define UINT_MAX __UINT_MAX__
#endif

/**
   @def GSL_RNG
   @brief gsl's mersenne twister rng
   @def KISS_RNG
   @brief Marsaglia's keep it simple stupid rng.
   @def MWC_1038_RNG
   @brief Marsaglia's ridiculous period generator
   @def MWC_4096_RNG
   @brief Marsaglia's ridiculous(er) period generator
   @def WELL_RNG
   @brief the default, well equidistributed generator. Improved MT

   @warning defaults to the WELL_RNG
 */

// Default to the WELL RNG if gsl is specified but no header found
#if ( defined GSL_RNG ) && ( defined HAVE_GSL )  
  #include <gsl/gsl_rng.h>
  #define RNG_TABLE 1
#elif defined KISS_RNG
  #define RNG_TABLE 4
#elif defined MWC_1038_RNG
  #define RNG_TABLE 1038
#elif defined MWC_4096_RNG
  #define RNG_TABLE 4096
#else
  #define WELL_RNG
  #define RNG_TABLE 624
#endif

/******************************

    Control for the U(1) code

*******************************/

/**
   @def U1_DFT
   @brief U(1) fourier transform code, defaults to the slower dft because I haven't satisfactorily proved the DHT is equivalent. I believe it is and it is about twice as fast
 */
#ifndef U1_DHT
  #define U1_DFT
#endif

#endif
