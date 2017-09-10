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

// computes the first 12 moments of Q in T-direction
static int
Time_Moments( const GLU_complex *qtop ,
	      const char *str ,
	      const double NORM )
{
  // sum qtop into a correlator q(t)
  GLU_complex *Qt = malloc( Latt.dims[ ND-1 ] * sizeof( GLU_complex ) ) ;
  
  // compute the moments \sum_t (-1)^n t^{2n}/(2n!) Q(t)
  const double fac[ 12 ] = { 1 , 1 , 2 , 6 , 24 , 120 , 720 ,
			     5040 , 40320 , 362880 , 362800 ,
			     39916800 } ;
  const int LT = (int)Latt.dims[ ND-1 ] ;
  int n , t ;

  // compute the temporal sums
  #pragma omp parallel for private(t)
  for( t = 0 ; t < Latt.dims[ ND-1 ] ; t++ ) {
    Qt[ t ] = 0.0 ;
    register GLU_complex sum = 0.0 ;
    size_t i ;
    for( i = t*LCU ; i < (t+1)*LCU ; i++ ) {
      sum += qtop[i] ;
    }
    Qt[t] = sum * NORM ;
  }

  // initialise powers of t to t^0
  int *tpow = malloc( Latt.dims[ ND-1 ] * sizeof( int ) ) ;
  for( t = 0 ; t < LT ; t++ ) {
    tpow[ t ] = 1 ;
  }
  
  // loop increasing powers of moments of Q
  for( n = 0 ; n < 12 ; n++ ) {
    GLU_complex sum = 0.0 ;
    for( t = -LT/2 ; t < LT/2 ; t++ ) {
      const int posit = ( t + LT ) % LT ;
      sum += tpow[ posit ] * Qt[ posit ] ;
      tpow[ posit ] *= t ;
    }
    // multiply by correct i-factor, shift of Z_4
    switch( n&3 ) {
    case 0 : sum *= +1. / fac[ n ] ; break ;
    case 1 : sum *= +I  / fac[ n ] ; break ;
    case 2 : sum *= -1. / fac[ n ] ; break ;
    case 3 : sum *= -I  / fac[ n ] ; break ;
    }
    // divide by the factorial and sign
    fprintf( stdout , "[QMOMENTS] %s Moment_%d %1.12e %1.12e \n" ,
	     str , n ,
	     creal( sum ) , cimag( sum ) ) ;
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
compute_Q_moments( GLU_complex *qtop )
{
  // normalisations
  const double NORM = -0.001583143494411527678811 ; // -1.0/(64*Pi*Pi)
  
  Time_Moments( qtop , "Q" , NORM ) ;

  return GLU_SUCCESS ;
}

// compute the moments of the topological charge squared
// fftw will use qtop in here!
static int
compute_Q2_moments( GLU_complex *qtop )
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

  Time_Moments( qtop , "Q^2" , mulfac ) ;

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

  Time_Moments( out , "Q^2" , NORMSQ ) ;

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
  if( compute_Q_moments( qtop ) == GLU_FAILURE ) {
    goto memfree ;
  }

  // compute moments of Q^2
  if( compute_Q2_moments( qtop ) == GLU_FAILURE ) {
    goto memfree ;
  }

  FLAG = GLU_SUCCESS ;
  
 memfree :
  free( qtop ) ;
  
  return FLAG ;
}
