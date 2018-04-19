/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (Qcorr.c) is part of GLU.

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
   @file Qcorr.c
   @brief computation of the topological susceptibility correlator
 */
#include "Mainfile.h"

#include "clover.h"       // computation of the topological charge
#include "cut_output.h"   // automatic formatting of our output file
#include "cut_routines.h" // momentum cuts for config-space vector
#include "geometry.h"     // for the spacing computation
#include "GLU_bswap.h"    // byte swapping
#include "plan_ffts.h"    // config space correlator is convolution
#include "SM_wrap.h"      // in case we wish to do smearing

// compute the correlator using the convolved topological 
// charge correlator
static int
temporal_correlator( const GLU_complex *qtop ,
		     const char *str )
{
  // append .tcorr onto the output
  char cstr[ 1024 ] ;
  sprintf( cstr , "%s.tcorr" , str ) ; 
  FILE *Ap = fopen( cstr , "wb" ) ; 

  // timeslice length
  size_t lt[ 1 ] = { Latt.dims[ND-1] } ;

  double *ct = malloc( Latt.dims[ ND-1 ] * sizeof( double ) ) ;

  // write out the timeslice list ...
  write_tslice_list( Ap , lt ) ;

  // compute ct
  size_t t ;
#pragma omp parallel for private(t)
  for( t = 0 ; t < Latt.dims[ ND-1 ] ; t++ ) {
    const GLU_complex *p = qtop + LCU * t ;
    register double sum = 0.0 ;
    size_t i ;
    for( i = 0 ; i < LCU ; i++ ) {
      sum += creal( *p ) , p++ ;
    }
    ct[ t ] = sum ;
  }

  for( t = 0 ; t < Latt.dims[ ND-1 ] ; t++ ) {
    printf( "[QSUSC] CORR %zu %e \n" , t , ct[t] ) ;
  }

  // tell us where to go
  fprintf( stdout , "[CUTS] Outputting correlator to %s \n" , cstr ) ;

  // and write the props ....
  write_g2_to_list( Ap , ct , lt ) ;

  // memory frees
  fclose( Ap ) ;
  free( ct ) ;

  return GLU_SUCCESS ;
}

// computation of the correlator
// C(r) = < q( x ) q( y ) >  -> r = ( x - y ), r^2 < CUTINFO.max_mom
// over all indices, q(x) is the topological charge at site x
int
compute_Qcorr( struct site *__restrict lat ,
	       const struct cut_info CUTINFO ,
	       const size_t measurement )
{
  // normalisations
  const double NORM = -0.001583143494411527678811 ; // -1.0/(64*Pi*Pi)
  const double NORMSQ = NORM * NORM ;

  // info
  check_psq( CUTINFO ) ;

  // compute the ratios of the dimensions in terms of the smallest
  simorb_ratios( ND ) ;

  // if we already have a file, read it
  size_t size[1] = { 0 } ;
  struct veclist *list = compute_veclist( size , CUTINFO , ND , GLU_TRUE ) ;

  // set up the outputs
  char *str = output_str_struct( CUTINFO ) ;
  sprintf( str , "%s.m%zu" , str , measurement ) ;
  FILE *Ap = fopen( str , "wb" ) ; 

  // flag for success or failure
  int FLAG = GLU_FAILURE ;

  // write the list, in cut_outputs.c
  write_mom_veclist( Ap , size , list , ND ) ;

  // set up the matrix-valued array of the topological charge
  struct fftw_small_stuff FFTW ;
#ifndef HAVE_FFTW3_H
  FFTW.in = malloc( LVOLUME * sizeof( GLU_complex ) ) ;
#endif
  size_t i ;

  // allocate the results
  double *qcorr = malloc( size[0] * sizeof( double ) ) ; 

  // logic to fail
  if( list == NULL ) goto memfree ;
  if( qcorr == NULL ) goto memfree ;
  if( FFTW.in == NULL ) goto memfree ;

  // precompute all of charge densities q(x)
  compute_Gmunu_array( FFTW.in , lat ) ;

  register double sum = 0.0 , sumsq = 0.0 ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    sum += creal( FFTW.in[i] ) ;
    sumsq += creal( FFTW.in[i]*FFTW.in[i] ) ;
  }
  fprintf( stdout , "\n[QTOP] Q %zu %1.12e %1.12e %1.12e \n" ,
	   measurement , sum * NORM , sum * sum * NORMSQ , sumsq * NORMSQ ) ;
  
  // if we have fftw, use it for the contractions
#ifdef HAVE_FFTW3_H

  const double mulfac = NORMSQ / (double)LVOLUME ;

  // FFT Gmunu
  small_create_plans_DFT( &FFTW , Latt.dims , ND ) ;

  // computes the full correlator in p-space and FFTs back
  fftw_execute( FFTW.forward ) ;

  // is a convolution, Volume norm is for the FFT
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    FFTW.out[ i ] *= conj( FFTW.out[ i ] ) * mulfac ;
  }

  fftw_execute( FFTW.backward ) ;

  #ifdef verbose
  fprintf( stdout , "\n[QSUSC] Check sum :: %f \n" , creal( FFTW.in[0] ) ) ;
  #endif

  // loop the possible rsqs
  #pragma omp parallel for private(i)
  for( i = 0 ; i < size[0] ; i++ ) {
    qcorr[i] = (double)creal( FFTW.in[ list[i].idx ] ) ;
  }

  // warning, this code is super slow compared to the FFT convolution one
  // above it probably shouldn't be used unless under duress
#else

  // loop the possible rsqs
  #pragma omp parallel for private(i)
  for( i = 0 ; i < size[0] ; i++ ) {
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
      sumqq += creal( FFTW.in[source] * FFTW.in[sink] ) ;
    }
    qcorr[i] = sumqq * NORMSQ ;
  }
#endif

  // look at the correlator of this thing too
  temporal_correlator( FFTW.in , str ) ;

  // tell us where to go
  fprintf( stdout , "[CUTS] Outputting to %s \n" , str ) ;

  // write out the result of all that work
  write_g2_to_list( Ap , qcorr , size ) ;
  
  // we are quite successful if we get here
  FLAG = GLU_SUCCESS ;

 memfree :

  // close the file and its name
  fclose( Ap ) ;
  free( str ) ;

  // free up all of the allocations
  free( qcorr ) ;

#ifdef HAVE_FFTW3_H
  small_clean_up_fftw( FFTW ) ;
#else
  free( FFTW.in ) ;
#endif
  
  // free the list
  free( list ) ;

  return FLAG ;
}
