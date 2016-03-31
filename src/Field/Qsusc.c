/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (Qsusc.c) is part of GLU.

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
   @file Qsusc.c
   @brief computation of the topological susceptibility correlator
   @warning calls the smearing wrapper, which overwrites the links
 */

#include "Mainfile.h"

#include "clover.h"       // computation of the topological charge
#include "cut_output.h"   // automatic formatting of our output file
#include "cut_routines.h" // momentum cuts for config-space vector
#include "geometry.h"     // for the spacing computation
#include "GLU_bswap.h"    // byte swapping
#include "plan_ffts.h"    // config space correlator is convolution
#include "SM_wrap.h"      // in case we wish to do smearing

// computation of the correlator
// C(r) = < q( x ) q( y ) >  -> r = ( x - y ), r^2 < CUTINFO.max_mom
// over all indices, q(x) is the topological charge at site x
int
compute_Qsusc( struct site *__restrict lat ,
	       const struct cut_info CUTINFO ,
	       const struct sm_info SMINFO )
{
  // normalisations
  const double NORM = -0.001583143494411527678811 ; // 1.0/(64*Pi*Pi)
  const double NORMSQ = NORM * NORM ;

  // do some smearing ...
  SM_wrap_struct( lat , SMINFO ) ;

  // info
  check_psq( CUTINFO ) ;

  // compute the ratios of the dimensions in terms of the smallest
  simorb_ratios( ND ) ;

  // if we already have a file, read it
  int size[1] = { 0 } ;
  struct veclist *list = compute_veclist( size , CUTINFO , ND , GLU_TRUE ) ;

  // set up the outputs
  const char *str = output_str_struct( CUTINFO ) ;  
  FILE *Ap = fopen( str , "wb" ) ; 

  // flag for success or failure
  int FLAG = GLU_FAILURE ;

  // write the list, in cut_outputs.c
  write_mom_veclist( Ap , size , list , ND ) ;

  // set up the matrix-valued array of the topological charge
  GLU_complex *qtop = malloc( LVOLUME * sizeof( GLU_complex ) ) ;
  size_t i ;

  // allocate the results
  double *qcorr = malloc( size[0] * sizeof( double ) ) ; 

#ifdef verbose
  register double sum = 0.0 , sumsq = 0.0 ;
#endif

  // fft'd list and plans
#ifdef HAVE_FFTW3_H
  GLU_complex *out = NULL ;
  fftw_plan forward , backward ;
#endif

  // logic to fail
  if( list == NULL ) goto memfree ;
  if( qcorr == NULL ) goto memfree ;
  if( qtop == NULL ) goto memfree ;

  // precompute all of charge densities q(x)
  compute_Gmunu_array( qtop , lat ) ;

#ifdef verbose
  for( i = 0 ; i < LVOLUME ; i++ ) {
    sum += creal( qtop[i] ) ;
    sumsq += creal( qtop[i] * qtop[i] ) ;
  }
  fprintf( stdout , "\nQTOP %f %f \n" , sum * NORM , sumsq * NORMSQ ) ;
#endif
  
  // if we have fftw, use it for the contractions
#ifdef HAVE_FFTW3_H

  // init parallel threads, maybe
  if( parallel_ffts( ) == GLU_FAILURE ) {
    fprintf( stderr , "[PAR] Problem with initialising the OpenMP "
	              "FFTW routines \n" ) ;
    goto memfree ;
  }

  // FFT Gmunu
  out = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;

  small_create_plans_DFT( &forward , &backward , qtop , out , ND ) ;
  fftw_execute( forward ) ;
  
  // computes the full correlator in p-space and FFTs back
  // is a convolution, Volume norm is for the FFT
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    out[ i ] *= conj( out[ i ] ) / (double)LVOLUME ;
  }

  fftw_execute( backward ) ;

  // cleanup and memory deallocate
  fftw_destroy_plan( backward ) ;
  fftw_destroy_plan( forward ) ;
  fftw_cleanup( ) ;
  #ifdef OMP_FFTW
  fftw_cleanup_threads( ) ;
  #endif
  free( out ) ;

  #ifdef verbose
  fprintf( stdout , "\n[QSUSC] Check sumsq :: %f \n" , 
	   creal( qtop[0] ) * NORMSQ ) ;
  #endif

  // loop the possible rsqs
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < size[0] ; i++ ) {
    qcorr[i] = (double)creal( qtop[ list[i].idx ] ) * NORMSQ ;
  }

  // warning, this code is super slow compared to the FFT convolution one
  // above it probably shouldn't be used unless under duress
#else

  // loop the possible rsqs
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < size[0] ; i++ ) {
    // some storage for the traces
    register double sumqq = 0.0 ;

    int separation[ ND ] ;
    get_mom_2piBZ( separation , list[i].idx , ND ) ;

    // loop the lattice varying source and sink with the correct separation
    size_t source = 0 ;
    for( source = 0 ; source < LVOLUME ; source++ ) {
      // translate the source from k by a vector separation
      const size_t sink = compute_spacing( separation , source , ND ) ;
      //trace of the products
      sumqq += creal( qtop[source] * qtop[sink] ) ;
    }
    qcorr[i] = sumqq * NORMSQ ;
  }

#endif

  // tell us where to go
  fprintf( stdout , "[CUTS] Outputting to %s \n" , str ) ;

  // write out the result of all that work
  write_g2_to_list( Ap , qcorr , size ) ;
  
  // we are quite successful if we get here
  FLAG = GLU_SUCCESS ;

 memfree :

  // close the file and its name
  fclose( Ap ) ;
  free( (void*)str ) ;

  // free up all of the allocations
  free( qcorr ) ;  
  free( qtop ) ;

  // free the list
  free( list ) ;

  return FLAG ;
}
