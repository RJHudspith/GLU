/*
    Copyright 2013-2016 Renwick James Hudspith

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
Amu_fields( struct site *__restrict A ,
	    const lie_field_def def )
{
  size_t i ;

  // callback for the log definition
  void (*log)( GLU_complex Q[ NCNC ] ,
	       const GLU_complex U[ NCNC ] ) = Hermitian_proj ;
  switch( def ) {
  case LINEAR_DEF :
    break ;
  case LOG_DEF : 
    log = exact_log_slow ; 
    break ;
  }

#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    struct spt_site temp ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      log( temp.O[mu] , A[i].O[mu] ) ;
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
  size_t i ;
#pragma omp parallel for private(i) reduction(+:sum)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex dA[ NCNC ] GLUalign ;
    zero_mat( dA ) ;
    size_t mu , elem , back ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      // can we not use navigation here?
      back = gen_shift( i , -mu - 1 ) ; 
      for( elem = 0 ; elem < NCNC ; elem++ ) {
	dA[ elem ] += A[i].O[mu][elem] - A[back].O[mu][elem] ;
      }
    }
    GLU_real tr ;
    trace_ab_herm( &tr , dA , dA ) ;
    sum = sum + (double)tr ;
  }
  fprintf( stdout , "[CUTS] GF accuracy :: %e \n" , sum / ( NC * LVOLUME ) ) ;
  return ;
}

#ifdef HAVE_FFTW3_H

// 
static int
spatial_correlator( struct site *__restrict A ,
		    double *gsp ,
		    const double gsnorm ,
		    const size_t spacing )
{
  if( spacing > Latt.dims[ND-1] ) { return GLU_FAILURE ; }
  // init parallel threads, maybe
  if( parallel_ffts( ) == GLU_FAILURE ) {
    fprintf( stderr , "[PAR] Problem with initialising "
	              "the OpenMP FFTW routines \n" ) ;
    // should clean up the memory here
    return GLU_FAILURE ;
  }

  // FFTW routines
  GLU_complex *out = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;
  GLU_complex *in = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;

  // create some plans
  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , in , out , ND ) ;

  //forward transform
  size_t mu , i , j ;
  for( mu = 0 ; mu < ND ; mu++ ) {
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
    #ifdef CUT_FORWARD
    out[i] = 0 ;
    #else
    in[i] = 0 ;
    #endif
    for( mu = 0 ; mu < ND ; mu++ ) {
      trace_ab_dag( &tr , A[i].O[mu] , A[i].O[mu] ) ;
      sum += tr ; 
    }
    #ifdef CUT_FORWARD
    out[i] = sum ;
    #else
    in[i] = sum ;
    #endif
  }

  // fft back into config space
#ifdef CUT_FORWARD
  fftw_execute( backward ) ;
#else
  fftw_execute( forward ) ;
#endif

  // compute the sum
  size_t t ;
#pragma omp parallel for private(i)
  for( t = 0 ; t < Latt.dims[ND-1] ; t++ ) {
    #ifdef CUT_FORWARD
    GLU_complex *p = in + LCU * t ;
    #else
    GLU_complex *p = out + LCU * t ;
    #endif
    register double sum = 0 ;
    for( i = 0 ; i < LCU ; i++ ) {
      sum += creal( *p ) ;
      p++ ;
    }
    gsp[ t ] = sum * gsnorm ;
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

  return GLU_SUCCESS ;
}

#else

// contract = Tr( A_\mu(x) A_\mu(y) ) ( for \mu < ND )
static double
contract_slices( const struct site *__restrict A ,
		 const size_t x ,
		 const size_t t_ref )
{
  GLU_real tr ;
  register double loc_tr = 0.0 ;
  size_t j ;
  for( j = 0 ; j < LCU ; j++ ) {
    const size_t y = j + t_ref ;
    #if ND == 4
    trace_ab_herm( &tr , A[x].O[0] , A[y].O[0] ) ;
    loc_tr += (double)tr ;
    trace_ab_herm( &tr , A[x].O[1] , A[y].O[1] ) ;
    loc_tr += (double)tr ;
    trace_ab_herm( &tr , A[x].O[2] , A[y].O[2] ) ;
    loc_tr += (double)tr ;
    //trace_ab_herm( &tr , A[x].O[3] , A[y].O[3] ) ;
    //loc_tr += (double)tr ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      trace_ab_herm( &tr , A[x].O[mu] , A[y].O[mu] ) ;
      loc_tr += (double)tr ;
    }
    #endif
  }
  return loc_tr / (double)( LCU ) ;
}

// recursion for the sum on a sparse wall type thing
static double
recurse_sum( const struct site *__restrict A ,
	     const size_t idx1 , 
	     const size_t idx2 ,
	     const size_t spacing ,
	     int vec[ ND ] , 
	     const size_t mu , 
	     size_t *norm ) 
{
  if( mu > ND-2 ) { 
    const size_t xidx = idx1 + gen_site( vec ) ;
    *norm = *norm + 1 ;
    return contract_slices( A , xidx , idx2 ) ;
  }

  double sum = 0.0 ;
  size_t nu ;
  int loc_vec[ ND ] ; 
  for( nu = 0 ; nu < ND ; nu++ ) {
    loc_vec[ nu ] = vec[ nu ] ;
  }
  for( nu = 0 ; nu < Latt.dims[mu] ; nu+=spacing ) {
    loc_vec[ mu ] = nu ;
    sum += recurse_sum( A , idx1 , idx2 , spacing , 
			loc_vec , mu + 1 , norm ) ;
  }
  return sum ;
}

/////////// spatial-spatial point sources over the whole lattice  ////////
static int
spatial_pointsource( const struct site *__restrict A ,
		     double *__restrict gs ,
		     const double gsnorm )
{
  // loop time 
  size_t t ;
#pragma omp parallel for private(t)
  for( t = 0 ; t < LT/2+1 ; t++ ) { // is exactly symmetric around LT/2
    // local trace
    register double loc_tr = 0.0 ;
    // temporary trace
    GLU_real tr ;
    // loop separations tau
    size_t tau , i , j ;
    for( tau = 0 ; tau < LT ; tau++ ) {
      const size_t idx1 = LCU * tau ;
      const size_t idx2 = LCU * ( ( tau + t ) % LT ) ;
      // loop spatial hypercube index i
      for( i = 0 ; i < LCU ; i++ ) {
	const size_t x = i + idx1 ;
	for( j = 0 ; j < LCU ; j++ ) {
	  const size_t y = j + idx2 ;
	  #if ND == 4
	  trace_ab_herm( &tr , A[x].O[0] , A[y].O[0] ) ;
	  loc_tr += (double)tr ;
	  trace_ab_herm( &tr , A[x].O[1] , A[y].O[1] ) ;
	  loc_tr += (double)tr ;
	  trace_ab_herm( &tr , A[x].O[2] , A[y].O[2] ) ;
	  loc_tr += (double)tr ;
	  //trace_ab_herm( &tr , A[x].O[3] , A[y].O[3] ) ;
	  //loc_tr += (double)tr ;
	  #else
	  size_t mu ;
	  for( mu = 0 ; mu < ND-1 ; mu++ ) {
	    trace_ab_herm( &tr , A[x].O[mu] , A[y].O[mu] ) ;
	    loc_tr += (double)tr ;
	  }
	  #endif
	}
      }
    }
    gs[t] = loc_tr * 2 / ( LCU * LCU * LT ) ;
  }
#pragma omp parallel for private(t)
  for( t = LT/2+1 ; t < LT ; t++ ) {
    // all to all is symmetric around LT/2
    gs[ t ] = gs[ (LT)-t ] ;
  }
  return GLU_SUCCESS ;
}

// compute the correlator from a specific point "x"
static int
spatial_correlator( const struct site *__restrict A ,
		    double *__restrict gs ,
		    const double gsnorm ,
		    const size_t spacing )
{
  // this one is more effective
  if( spacing == 1 ) {
    return spatial_pointsource( A , gs , gsnorm ) ;
  }
  // this one is more general
  size_t t ;
#pragma omp parallel for private(t)
  for( t = 0 ; t < LT ; t++ ) { 
    // local trace
    double loc_tr = 0.0 ;
    size_t norm = 0 , tau , mu ;
    for( tau = 0 ; tau < LT ; tau++ ) {
      const size_t idx1 = LCU * tau ;
      const size_t idx2 = LCU * ( ( tau + t ) % LT ) ;
      int x[ ND ] ;
      for( mu = 0 ; mu < ND ; mu++ ) { x[mu] = 0 ; }
      loc_tr += recurse_sum( A , idx1 , idx2 , spacing , 
			     x , 0 , &norm ) ;
    }
    gs[t] = loc_tr * 2.0 / (double)norm ;
  }
  // we are so successful
  return GLU_SUCCESS ;
}

#endif

// computes the configuration space propagators
int 
cuts_struct_configspace( struct site *__restrict A ,
			 const struct cut_info CUTINFO ,
			 const struct sm_info SMINFO )
{
  fprintf( stdout , "\n[CUTS] Computing the configuration-space"
	   "gluon propagator\n" ) ;

  // smeared gluon field to extract the ground state?
  SM_wrap_struct( A , SMINFO ) ;

  // take the log
  if( Amu_fields( A , CUTINFO.definition ) != GLU_SUCCESS ) {
    fprintf( stderr , "[CUTS] Something wrong in Logging of "
	              "the fields, check it out \n" ) ;
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

  // flag whether it worked or not
  int FLAG = GLU_FAILURE ;

  // timeslice length
  size_t lt[ 1 ] = { LT } ;

  // spatial correlator from the convolution in momentum space
  if( spatial_correlator( A , gsp , spat_norm , 
			  CUTINFO.max_t ) == GLU_FAILURE ) {
    goto memfree ;
  }

  // write out the timeslice list ...
  write_tslice_list( Ap , lt ) ;

  // and write the props ....
  write_g2_to_list( Ap , gsp , lt ) ;

  // tell us where we are writing out to
  fprintf( stdout , "[CUTS] data written to %s \n" , str ) ;

  // we are so successful
  FLAG = GLU_SUCCESS ;

 memfree :

  // close the file
  fclose( Ap ) ;
  // and its name
  free( str ) ;

  // and free the gluon correlator
  free( gsp ) ;

  return FLAG ;
}

#undef LT // clean it up?
