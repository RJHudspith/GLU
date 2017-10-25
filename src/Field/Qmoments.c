/*
    Copyright 2013-2017 Renwick James Hudspith

    This file (Qmoments.c) is part of GLU.

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
   @file Qmoments.c
   @brief compute the moments of the topological charge
 */
#include "Mainfile.h"

#include "clover.h"       // computation of the topological charge
#include "cut_output.h"   // automatic formatting of our output file
#include "cut_routines.h" // momentum cuts for config-space vector
#include "geometry.h"     // for the spacing computation
#include "GLU_bswap.h"    // byte swapping
#include "plan_ffts.h"    // config space correlator is convolution
#include "SM_wrap.h"      // in case we wish to do smearing

// number of moments that we look at
#define NQMOMENTS ((size_t)12)

// local version of the one in geometry
static void
get_mom_2piBZ_loc( int x[ ND ] ,
		   const size_t Dims[ ND ] ,
		   const size_t i ,
		   const size_t DIMS )
{
  int mu , subvol = 1 , sum = 0 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( mu != (int)DIMS ) {
      x[ mu ] = ( ( i - i % subvol ) / subvol ) % Dims[ mu ] ;
      subvol *= Dims[ mu ] ;
      sum += x[mu] ;
    } else {// set it to 0?
      x[ mu ] = 0 ;
    }
  }
  return ;
}

// given a local SLAB site index return the global index taking
// into account for the translational invariance term T0
static size_t
loc_idx( const size_t i ,
	 const size_t dims[ ND ] ,
	 const size_t DIR ,
	 const size_t slice )
{
  int x[ ND ] ; 
  get_mom_2piBZ_loc( x , dims , i , ND ) ;
  x[ DIR ] = slice ;

  #ifdef verbose
  printf( "%d %d %d %d \n" , x[0] , x[1] , x[2] , x[3] ) ;
  #endif
  
  return gen_site( x ) ;
}

// computes the first 12 moments of Q in DIR-direction \sum_t (-1)^n t^{2n}/(2n!) Q(t)
static int
Time_Moments( double *Moment ,
	      const GLU_complex *qtop ,
	      const double NORM ,
	      const size_t DIR ,
	      const size_t measurement )
{
  // sum qtop into a correlator q(t)
  GLU_complex *Qt = malloc( Latt.dims[ DIR ] * sizeof( GLU_complex ) ) ;
  
  // factorials
  const double fac[ NQMOMENTS ] = { 1 , 1 , 2 , 6 , 24 , 120 , 720 ,
				    5040 , 40320 , 362880 , 362800 ,
				    39916800 } ;
  const int LT = (int)Latt.dims[ DIR ] ;
  size_t dims[ ND ] ;
  int n , t , subvol = 1 ;

  for( n = 0 ; n < ND ; n++ ) {
    if( n != DIR ) {
      dims[ n ] = Latt.dims[ n ] ;
      subvol *= Latt.dims[ n ] ;
    } else {
      dims[ n ] = 1 ;
    }
  }

  // compute the temporal sums
#pragma omp parallel for private(t)
  for( t = 0 ; t < Latt.dims[ DIR ] ; t++ ) {
    Qt[ t ] = 0.0 ;
    register GLU_complex sum = 0.0 ;
    size_t i ;
    for( i = 0 ; i < subvol ; i++ ) {
      sum += qtop[ loc_idx( i , dims , DIR , t ) ] ;
    }
    Qt[t] = sum * NORM ;
  }

  // initialise powers of t to t^0
  int *tpow = malloc( Latt.dims[ DIR ] * sizeof( int ) ) ;
  for( t = 0 ; t < LT ; t++ ) {
    tpow[ t ] = 1 ;
  }
  
  // loop increasing powers of moments of Q
  for( n = 0 ; n < NQMOMENTS ; n++ ) {
    double sum = 0.0 ;
    for( t = -LT/2+1 ; t < LT/2 ; t++ ) {
      const int posit = ( t + LT ) % LT ;
      sum += tpow[ posit ] * creal( Qt[ posit ] ) ;
      tpow[ posit ] *= t ;
    }
    sum /= fac[ n ] ;
    
    // multiply by correct i-factor, shift of Z_4
    switch( n&3 ) {
    case 0 : Moment[ n + NQMOMENTS * DIR ] = sum  ; break ;
    case 1 : Moment[ n + NQMOMENTS * DIR ] = sum  ; break ;
    case 2 : Moment[ n + NQMOMENTS * DIR ] = -sum ; break ;
    case 3 : Moment[ n + NQMOMENTS * DIR ] = -sum ; break ;
    }
    
    // divide by the factorial and sign
    //#ifdef verbose
    fprintf( stdout , "[QMOMENTS] Dir_%zu Moment_%d %1.12e \n" ,
	     DIR , n , Moment[ n + NQMOMENTS * DIR ] ) ;
    //#endif
  }
  
  if( tpow != NULL ) {
    free( tpow ) ;
  }

  if( Qt != NULL ) {
    free( Qt ) ;
  }

  return GLU_SUCCESS ;
}

// compute the moments of the topological charge
static int
compute_Q_moments( struct Qmoments *Qmom ,
		   GLU_complex *qtop ,
		   const size_t measurement )
{
  // normalisations
  const double NORM = -0.001583143494411527678811 ; // -1.0/(64*Pi*Pi)

    // init parallel threads, maybe
  if( parallel_ffts( ) == GLU_FAILURE ) {
    fprintf( stderr , "[PAR] Problem with initialising the OpenMP "
	              "FFTW routines \n" ) ;
    return GLU_FAILURE ;
  }

  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Time_Moments( Qmom -> Q , qtop ,
		  NORM , mu , measurement ) ;
  }

  return GLU_SUCCESS ;
}

// compute the moments of the topological charge squared
// fftw will use qtop in here!
static int
compute_Q2_moments( struct Qmoments *Qmom ,
		    GLU_complex *qtop ,
		    const size_t measurement )
{
  GLU_complex *out = NULL ;
  const double NORM = -0.001583143494411527678811 ; // -1.0/(64*Pi*Pi)
  const double NORMSQ = NORM * NORM ;
  const double mulfac = NORMSQ / (double)LVOLUME ;
  size_t i ;
  
#ifdef HAVE_FFTW3_H
  // do the convolution of below
  out = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;

    // init parallel threads, maybe
  if( parallel_ffts( ) == GLU_FAILURE ) {
    fprintf( stderr , "[PAR] Problem with initialising the OpenMP "
	              "FFTW routines \n" ) ;
    return GLU_FAILURE ;
  }

  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , qtop , out , Latt.dims , ND ) ;

  // computes the full correlator in p-space and FFTs back
  fftw_execute( forward ) ;

  // is a convolution, Volume norm is for the FFT
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    out[ i ] *= conj( out[ i ] ) ;
  }

  fftw_execute( backward ) ;

  // cleanup and memory deallocate
  fftw_destroy_plan( backward ) ;
  fftw_destroy_plan( forward ) ;
  fftw_cleanup( ) ;
  #ifdef OMP_FFTW
  fftw_cleanup_threads( ) ;
  #endif

  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Time_Moments( Qmom -> Q2 , qtop ,
		  mulfac , mu , measurement ) ;
  }
  
#else

  out = malloc( LVOLUME * sizeof( GLU_complex ) ) ;
  
  // compute all shifts
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register GLU_complex sum = 0.0 ;
    int separation[ ND ] ;
    get_mom_2piBZ( separation , i , ND ) ;
    size_t source ;
    for( source = 0 ; source < LVOLUME ; source++ ) {
      const size_t sink = compute_spacing( separation , source , ND ) ;
      sum += ( qtop[source] * qtop[sink] ) ;
    }
    out[i] = sum ;
  }

  size_t ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Time_Moments( Qmom -> Q2 , out ,
		  NORMSQ , mu , measurement ) ;
  }
  
#endif

  // free temporary "out" matrix
  if( out != NULL ) {
    free( out ) ;
  }
  
  return GLU_SUCCESS ;
}

// compute the topological moments
int
compute_Qmoments( struct site *__restrict lat ,
		  struct Qmoments *Qmom ,
		  const struct cut_info CUTINFO ,
		  const size_t measurement )
{
  // normalisations
  const double NORM = -0.001583143494411527678811 ; // -1.0/(64*Pi*Pi)
  
  GLU_complex *qtop = malloc( LVOLUME * sizeof( GLU_complex ) ) ;
  int FLAG = GLU_FAILURE ;

  // precompute all of the charge densities q(x)
  compute_Gmunu_array( qtop , lat ) ;

  // do a check
  double sum = 0.0 ;
  size_t i , n ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    sum += creal( qtop[i] ) ;
  }
  for( n = 1 ; n < 7 ; n++ ) {
    printf( "[QMOMENTS] Q^%zu %1.12e \n" ,
	    n , pow( sum*NORM , n ) ) ;
  }

  // compute moments of Q
  if( compute_Q_moments( Qmom , qtop , measurement ) == GLU_FAILURE ) {
    goto memfree ;
  }

  // compute moments of Q^2
  if( compute_Q2_moments( Qmom , qtop , measurement ) == GLU_FAILURE ) {
    goto memfree ;
  }

  FLAG = GLU_SUCCESS ;
  
 memfree :
  free( qtop ) ;
  
  return FLAG ;
}
