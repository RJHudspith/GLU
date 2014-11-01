/*
    Copyright 2013 Renwick James Hudspith

    This file (config_gluons.c) is part of GLU.

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
   @file config_gluons.c
   @brief code to compute the configuration space gluonic correlation functions
 */

#include "Mainfile.h"

#include "cut_output.h"   // output file
#include "cut_routines.h" // momentum cutting?
#include "geometry.h"     // general lexi site
#include "plan_ffts.h"    // config space correlator is convolution
#include "SM_wrap.h"      // do we want some smearing?

/**
   @param LT
   @brief length of the temporal direction
 */
#define LT Latt.dims[ND-1]

// compute the log fields, overwrite the link matrices!
static int
Amu_fields( A , def )
     struct site *__restrict A ;
     const lie_field_def def ;
{
  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    struct spt_site temp ;
    int mu ;
    switch( def )
      {
      case LINEAR_DEF :
	for( mu = 0 ; mu < ND ; mu++ ) 
	  Hermitian_proj( temp.O[mu] , A[i].O[mu] ) ;
	break ;
      case LOG_DEF :
	for( mu = 0 ; mu < ND ; mu++ ) 
	  exact_log_slow( temp.O[mu] , A[i].O[mu] ) ;
	break ;
      }
    memcpy( &A[i] , &temp , sizeof( struct spt_site ) ) ;
  }
  return GLU_SUCCESS ;
}

// compute the gauge fixing accuracy
static void
check_landau_condition( const struct site *__restrict A )
{
  // compute the gauge fixing accuracy Tr( dA dA )
  double sum = 0.0 ;
  int i ;
#pragma omp parallel for private(i) reduction(+:sum)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex dA[ NCNC ] = {} ;
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      const int back = gen_shift( i , -mu - 1 ) ; 
      int elem ;
      for( elem = 0 ; elem < NCNC ; elem++ ) {
	dA[ elem ] += A[i].O[mu][elem] - A[back].O[mu][elem] ;
      }
    }
    GLU_real tr ;
    trace_ab_herm( &tr , dA , dA ) ;
    sum = sum + (double)tr ;
  }
  printf( "[CUTS] GF accuracy :: %e \n" , 4.0 * sum / ( NC * LVOLUME ) ) ;
  return ;
}

// 
static void
spatial_corr_convolve( struct site *__restrict A ,
		       double *gsp ,
		       const double spat_norm )
{
  // FFTW routines
  GLU_complex *out = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;
  GLU_complex *in = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;

  // create some plans
  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , in , out , ND ) ;

  //forward transform
  int mu , i ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    int j ;
    for( j = 0 ; j < NCNC ; j++ ) {
      // FORWARD ONE
      #ifdef CUT_FORWARD
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	in[i] = A[i].O[mu][j] ;
      }
      fftw_execute( forward ) ;
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	A[i].O[mu][j] = out[i] ;
      }
      #else
      // BACKWARD TRANSFORM
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	out[i] = A[i].O[mu][j] ;
      }
      fftw_execute( backward ) ;
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	A[i].O[mu][j] = in[i] ;
      }
      #endif
    }
  }

  // compute the convolution in momentum space
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex tr , sum = 0.0 ;
    out[i] = 0 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      trace_ab_dag( &tr , A[i].O[mu] , A[i].O[mu] ) ;
      sum += tr ; 
    }
    #ifdef CUT_FORWARD
    out[i] = out[i] + sum ;
    #else
    in[i] = in[i] + sum ;
    #endif
  }

  // fft back into config space
#ifdef CUT_FORWARD
  fftw_execute( backward ) ;
#else
  fftw_execute( forward ) ;
#endif

  // compute the sum
  int t ;
#pragma omp parallel for private(i)
  for( t = 0 ; t < LT ; t++ ) {
    #ifdef CUT_FORWARD
    GLU_complex *p = in + LCU * t ;
    #else
    GLU_complex *p = out + LCU * t ;
    #endif
    register double sum = 0 ;
    for( i = 0 ; i < LCU ; i++ ) {
      sum += *p ;
      p++ ;
    }
    gsp[ t ] = sum * spat_norm ;
    printf( "%d %f \n" , t , gsp[t] ) ;
  }

  // free the FFTs
  fftw_destroy_plan( backward ) ;
  fftw_destroy_plan( forward ) ;
  fftw_cleanup( ) ;
#ifdef OMP_FFTW
  fftw_cleanup_threads( ) ;
#endif
  fftw_free( out ) ;  
  fftw_free( in ) ;
  return ;
}

// computes the configuration space propagators
int 
cuts_struct_configspace( struct site *__restrict A ,
			 const struct cut_info CUTINFO ,
			 const struct sm_info SMINFO )
{
  printf( "\n[CUTS] Computing the configuration-space gluon propagator\n" ) ;

  // smeared gluon field to extract the ground state?
  SM_wrap_struct( A , SMINFO ) ;

  // take the log
  if( Amu_fields( A , CUTINFO.definition ) != GLU_SUCCESS ) {
    printf( "[CUTS] Something wrong in Logging of the fields, check it out \n" ) ;
    return GLU_FAILURE ;
  }

  // check that we are in Landau gauge
  check_landau_condition( A ) ;

  // set up the outputs
  char *str = output_str_struct( CUTINFO ) ;  
  FILE *Ap = fopen( str , "wb" ) ; 

  // compute the temporal and spatial correlators
  double *gsp = malloc( Latt.dims[ ND-1 ] * sizeof( double ) ) ;

  // normalisations :: LVOLUME is the FFT norm
  const double spat_norm = 1.0 / ( ( ND-1 ) * LCU * LVOLUME ) ;
  
  // spatial correlator from the convolution in momentum space
  spatial_corr_convolve( A , gsp , spat_norm ) ;

  // write out the timeslice list ...
  int lt[ 1 ] = { LT } ;
  write_tslice_list( Ap , lt ) ;

  // and write the props ....
  write_g2_to_list( Ap , gsp , lt ) ;

  // tell us where we are writing out to
  printf( "[CUTS] data written to %s \n" , str ) ;

  // and free the gluon correlator
  free( gsp ) ;

  return GLU_SUCCESS ;
}

#undef LT // clean it up?
