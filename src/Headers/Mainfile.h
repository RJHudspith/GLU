/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (Mainfile.h) is part of GLU.

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
   @file Mainfile.h
   @brief the header that most files include, has all the definitions and matrix op libraries
 */
#ifndef GLU_MAINFILE_H
#define GLU_MAINFILE_H

// firstly this include
#include <config.h>

// are we compiling with either gsl or lapacke 
// only needed for certain SU(Nc>3) routines
#if (defined HAVE_GSL && defined HAVE_LAPACKE_H )
  #undef HAVE_GSL // lapacke routines are much faster!
#endif

// stop having a mixture of omp and cilk routines default
// to openmp routines because of better documentation
#if (defined _OPENMP ) && (defined HAVE_OMP_H ) && (defined HAVE_CILK_H)
  #undef HAVE_CILK_H
#endif

// wrap these openmp functions
#if (defined _OPENMP ) && (defined HAVE_OMP_H )
  #include <omp.h>
  #define get_GLU_thread() omp_get_thread_num()
#else
  #define get_GLU_thread() (0)
#endif

// generic includes ...
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>

// if we have FFTW3 we include it everywhere
#ifdef HAVE_FFTW3_H
#include <fftw3.h>
#endif

// included for the standard defintions of uint64_t, uint32_t and uint16_t
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

/**
   @def likely( x )
   @brief a gcc builtin expect tries harder in ifs when this is triggered

   @def unlikely( x )
   @brief a gcc builtin expect, not expecting the if to return true in this case
 */
// some gcc-oonly macros should guard these?
// #if ( __GNUC__ > 3 ) or something ? Well, xlc for the q and icc gets these
  #define likely( x )    __builtin_expect( ( x ) , 1 )
  #define unlikely( x )  __builtin_expect( ( x ) , 0 )

// include the definitions
#include "GLU_definitions.h"

// include the enumerated types ...
#include "GLU_enum.h"

// include the struct definitions that I use
#include "GLU_types.h"

/**
   @var Latt
   @brief extern the Lattice information "Latt" to everyone
   this struct is defined in GLU_types.h
 */
extern struct latt_info Latt ;

// my own necessary libs ...
#include "U_Nops.h"                // many matrix operations
#include "trace_abc.h"             // trace of the product of 3 matrices
// matrix multiplication
#include "MMUL.h"                  // a = b x c :: b,c
#include "MMUL_SUNC.h"             // a = b x c :: b,c in SU(NC)
#include "MMUL_dag.h"              // a = b x c^{\dagger}
#include "MMUL_dag_SUNC.h"         // a = b x c^{\dagger} SU(NC) variant
#include "MMULdag.h"               // a = b^{\dagger} x c 
#include "MMULdag_SUNC.h"          // a = b^{\dagger} x c SU(NC) variant
#include "MMULdagdag.h"            // a = b^{\dagger} x c^{\dagger}
#include "MMULdagdag_SUNC.h"       // a = b^{\dagger} x c^{\dagger}
// logs and hermitian projections
#include "exactQ.h"                // Q = log( U ) (approximations+exact)
#include "expMat.h"                // U = exp( Q ) (exact)
#include "GLU_malloc.h"            // GLUey allocations

#endif
