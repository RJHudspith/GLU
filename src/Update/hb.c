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

// just update the diagonal
//#define DIAGONAL_UPDATE

// hit the subgroups at random
//#define STOCH (NC*NC/4)

// maximum number of heatbath updates
#define NHBMAX (42)

// Creutz heatbath algorithm
static int
Creutz( double *__restrict d ,
	const double xl ,
	const double NORM ,
	const uint32_t thread ) 
{
  const register double x4 = par_rng_dbl( thread ) ;
  register const double a0 = 1.0 + \
    log( xl + ( 1.0 - xl ) * par_rng_dbl( thread ) ) * NORM ;
  *d = 1.0 - a0 ;
  return ( 1.0 - a0*a0 ) < x4*x4 ;
}

// Kennedy-Pendleton heatbath algorithm
static int
KP( double *__restrict d ,
    const double xl ,
    const double NORM ,
    const uint32_t thread ) 
{  
  const register double x4 = par_rng_dbl( thread ) ;
  const register double x3 = cos( TWOPI * par_rng_dbl( thread ) ) ;
  *d = -( log( par_rng_dbl( thread ) ) + log( par_rng_dbl( thread ) ) * x3*x3 ) * NORM ;
  return ( 1. - 0.5*(*d) ) < x4*x4 ;
}

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

  // for small beta (large NORM) the Creutz algorithm is preferred
  char msg[8] ;
  if( NORM < 2.0 ) {
    K = KP ; sprintf( msg , "KenPen" ) ;
  } else {
    K = Creutz ;
    xl = exp( -2./NORM ) ; 
    sprintf( msg , "Creutz" ) ;
  }

  // compute d which is (1-a0) = the leftmost su2 matrix element
  if( K( &d , xl , NORM , thread ) ) {
    // redo the algorithm
    while( K( &d , xl , NORM , thread ) && iters < NHBMAX ) {
      iters++ ;
    }
    // complain if we do too many hits?
    if( iters == NHBMAX ) {
      fprintf( stderr , "[HB] %d iterations of %s to no avail \n" ,
	       NHBMAX , msg ) ;
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

#ifdef DIAGONAL_UPDATE
  for( i = 0 ; i < NC-1 ; i++ ) {
    only_subgroup( &s0 , &s1 , &scale , U , staple , i ) ;
    if( generate_SU2( &s0 , &s1 , invbeta*scale , thread ) == GLU_FAILURE ) {
      continue ;
    }
    su2_rotate( U , s0 , s1 , i ) ;
  }
#elif defined STOCH
  for( i = 0 ; i < STOCH ; i++ ) {
    const size_t stoch = (size_t)( par_rng_dbl( thread ) * NSU2SUBGROUPS ) ;
    only_subgroup( &s0 , &s1 , &scale , U , staple , stoch ) ;
    if( generate_SU2( &s0 , &s1 , invbeta*scale , thread ) == GLU_FAILURE ) {
      continue ;
    }
    su2_rotate( U , s0 , s1 , stoch ) ;
  }
#else
  for( i = 0 ; i < NSU2SUBGROUPS ; i++ ) {
    only_subgroup( &s0 , &s1 , &scale , U , staple , i ) ;
    if( generate_SU2( &s0 , &s1 , invbeta*scale , thread ) == GLU_FAILURE ) {
      continue ;
    }
    su2_rotate( U , s0 , s1 , i ) ;
  }
#endif
  return ;
}

// perform a heat-bath over the whole lattice
int
hb_lattice( struct site *lat ,
	    const double invbeta ,
	    const struct draughtboard db )
{
  size_t mu , i ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // single node until I figure out the coloring
#ifdef IMPROVED_SMEARING
#pragma omp single
    {
      for( i = 0 ; i < LVOLUME ; i++ ) {
	GLU_complex stap[ NCNC ] GLUalign ;
	zero_mat( stap ) ;
	all_staples_improve( stap , lat , i , mu , ND , SM_APE ) ;
	hb( lat[ i ].O[mu] , stap , invbeta , get_GLU_thread() ) ;
      }
    }
#else
    size_t c ;
    // loop draughtboard coloring
    for( c = 0 ; c < db.Ncolors ; c++ ) {
      // parallel loop over all sites with this coloring
#pragma omp for private(i) //schedule(dynamic)
      for( i = 0 ; i < db.Nsquare[c] ; i++ ) {
	const size_t it = db.square[c][i] ;
	GLU_complex stap[ NCNC ] GLUalign ;
	zero_mat( stap ) ;
	all_staples( stap , lat , it , mu , ND , SM_APE ) ;
	hb( lat[ it ].O[mu] , stap , invbeta , get_GLU_thread() ) ;
      }
      // and that is it
    }
  }
#endif
  return GLU_SUCCESS ;
}
