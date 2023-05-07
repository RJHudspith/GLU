/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (hb.c) is part of GLU.

    GLU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

a    GLU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GLU.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   @file hb.c
   @brief heatbath code

   Notes :: s0 and s1 are the top row of the SU(2) matrix, this is all we need
            due to the funny nature of links we have to loop over directions and 
	    then volume, this is a bit out of order for our AoS setup but fuck it

   Credit :: Lots of credit to T deGrand whose code (in MILC) I used as the basis 
             for this work and whose code I tested this against.
 */
#define _GNU_SOURCE // sincos

#include "Mainfile.h"

#include "par_rng.h"    // parallel rng
#include "staples.h"    // all_staples()
#include "SU2_rotate.h" // rotation
#include "relax.h"      // OR

#include "gramschmidt.h"
#include "GLU_malloc.h"
#include "geometry.h"

static struct site *lat2 = NULL ;

// maximum number of heatbath updates
#define NHBMAX (42)

// Tony and Brian's algorithm is better - sorry LatticeGuy!
#define USE_KP

#ifdef USE_KP
// Kennedy-Pendleton heatbath algorithm
static int
KP( double *__restrict d ,
    const double xl ,
    const double NORM ,
    const uint32_t thread ) 
{  
  register const double x4 = par_rng_dbl( thread ) ;
  register const double x3 = cos( TWOPI * par_rng_dbl( thread ) ) ;
  *d = -( log( par_rng_dbl( thread ) ) + log( par_rng_dbl( thread ) ) * x3*x3 ) * NORM ;
  return ( 1. - 0.5*(*d) ) < x4*x4 ;
}
#else
// Creutz heatbath algorithm
static int
Creutz( double *__restrict d ,
	const double xl ,
	const double NORM ,
	const uint32_t thread ) 
{
  register const double x4 = par_rng_dbl( thread ) ;
  register const double a0 = 1.0 + \
    log( xl + ( 1.0 - xl ) * par_rng_dbl( thread ) ) * NORM ;
  *d = 1.0 - a0 ;
  return ( 1.0 - a0*a0 ) < x4*x4 ;
}
#endif

// generate SU(2) matrix proportional to boltzmann weight
static int
generate_SU2( GLU_complex *s0 ,
	      GLU_complex *s1 ,
	      const double NORM ,
	      const int thread )
{
  double xl = 1.0 , d = 0.0 ;
  size_t iters = 1 ;

  // function pointer for the algorithm we want
  static int (*K)( double *__restrict d ,
		   const double xl ,
		   const double NORM ,
		   const uint32_t thread ) ;
#ifdef USE_KP
  K = KP ;
#else
  K = Creutz ;
  xl = exp( -2./NORM ) ; 
#endif

  // compute d which is (1-a0) = the leftmost su2 matrix element
  if( K( &d , xl , NORM , thread ) ) {
    // redo the algorithm
    while( K( &d , xl , NORM , thread ) && iters < NHBMAX ) {
      iters++ ;
    }
  }

  // generate SU(2) matrix
  const double a0  = 1.0 - d ;
  const double ar2 = fabs( 1.0 - a0*a0 ) ; 
  const double a3  = sqrt( ar2 ) * ( 2.0 * par_rng_dbl( thread ) - 1.0 ) ;
  const double rho = sqrt( fabs( ar2 - a3*a3 ) ) ;
  const double xr2 = TWOPI * par_rng_dbl( thread ) ;

  // use sincos as it is cheaper
  double a1 , a2 ;
  sincos( xr2 , &a2 , &a1 ) ;
  a1 *= rho ; 
  a2 *= rho ;

  const double rS0 = creal( *s0 ) , iS0 = cimag( *s0 ) ;
  const double rS1 = creal( *s1 ) , iS1 = cimag( *s1 ) ;

  // unrolled matrix multiply
  *s0 = a0 * rS0 + a3 * iS0 + a2 * rS1 + a1 * iS1 + 
    I * ( -a0 * iS0 + rS0 * a3 - a2 * iS1 + rS1 * a1 ) ;
  *s1 = -a0 * rS1 + a3 * iS1 + a2 * rS0 - a1 * iS0 +
    I * ( -a0 * iS1 - a3 * rS1 + a2 * iS0 + a1 * rS0 ) ;

  return GLU_SUCCESS ;
}

// update matrices compute SU(2) subgroup of the staple
// compute heat-bath update and multiply through link
static void
hb( GLU_complex U[ NCNC ] , 
    const GLU_complex staple[ NCNC ] ,
    const double invbeta ,
    const uint32_t thread )
{
  GLU_complex s0 GLUalign , s1 GLUalign ;
  double scale GLUalign ;
  size_t i ;
  for( i = 0 ; i < NSTOCH ; i++ ) {
    #if (NSTOCH != NSU2SUBGROUPS)
    const size_t stoch = (size_t)( par_rng_dbl( thread ) * NSU2SUBGROUPS ) ;
    #else
    const size_t stoch = i ;
    #endif
    only_subgroup( &s0 , &s1 , &scale , U , staple , stoch ) ;
    generate_SU2( &s0 , &s1 , invbeta*scale , thread ) ;
    su2_rotate( U , s0 , s1 , stoch ) ;
  }
  return ;
}

// perform a heat-bath over the whole lattice
int
hb_lattice( struct site *lat ,
	    const double invbeta ,
	    const struct draughtboard db )
{
  // single node until I figure out the coloring -- could this just be a simple copy? i.e. a temporary gauge field being set and then pointer-swapped?
#ifdef IMPROVED_SMEARING
  size_t cmu , i ;
  for( cmu = 0 ; cmu < db.Ncolors ; cmu++ ) {
    const size_t mu = (size_t)cmu/4 ;
    // parallel loop over all sites with this coloring
    #pragma omp for private(i) 
    for( i = 0 ; i < db.Nsquare[cmu] ; i++ ) {
      const size_t it = db.square[cmu][i] ;
      GLU_complex stap[ NCNC ] GLUalign ;
      zero_mat( stap ) ;
      all_staples_improve( stap , lat , it , mu , ND , SM_APE ) ;
      hb( lat[ it ].O[mu] , stap , invbeta , get_GLU_thread() ) ;
      gram_reunit( lat[it].O[mu] ) ;
    }
  }
#else
  size_t cmu , i ;
  for( cmu = 0 ; cmu < ND*db.Ncolors ; cmu++ ) {
    const size_t mu = cmu/db.Ncolors ;
    const size_t c  = cmu%db.Ncolors ;
    // parallel loop over all sites with this coloring
    #pragma omp for private(i) 
    for( i = 0 ; i < db.Nsquare[c] ; i++ ) {
      const size_t it = db.square[c][i] ;
      GLU_complex stap[ NCNC ] GLUalign ;
      zero_mat( stap ) ;
      all_staples( stap , lat , it , mu , ND , SM_APE ) ;
      hb( lat[ it ].O[mu] , stap , invbeta , get_GLU_thread() ) ;
      gram_reunit( lat[it].O[mu] ) ;
    }
  }
#endif
  return GLU_SUCCESS ;
}

#undef NHBMAX
