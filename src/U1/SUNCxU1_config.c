/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (su3xu1_config.c) is part of GLU.

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
   @file SUNCxU1_config.c
   @brief U(1)-ifies the link matrices

   @warning requires FFTW linking to work
 */

#include "Mainfile.h"

// and the other headers it uses
#include "geometry.h"
#include "par_rng.h"
#include "GLU_timer.h"
#include "plan_ffts.h"
#include "U1_obs.h"

//#define U1_DHT

// just to make it clear what we are doing
enum{ CONJUGATE_NOT_IN_LIST , CONJUGATE_IN_LIST } ;

#ifdef HAVE_FFTW3_H

// little inline for the ( 0 , 0 , .. , 0 ) point in the -Pi -> Pi BZ
static size_t
compute_zeropoint( void )
{
#if ND == 4
  return ( Latt.dims[0] * ( 1 + Latt.dims[1] * ( 1 + Latt.dims[2] * ( 1 + Latt.dims[3] ) ) ) ) >> 1 ;
#else
  size_t check = Latt.dims[ ND-1 ] , mu ;
  for( mu = ND-2 ; mu != 0 ; mu-- ) {
    check = Latt.dims[ mu ] * ( 1 + check ) ;
  }
  return check >> 1 ;
#endif
}

// compute the conjugate site to "i"
static size_t
conjugate_site( const size_t i )
{
  int x[ ND ] ;
  get_mom_pipi( x , i , ND ) ;
  #if ND == 4
  x[0] = -x[0] ;
  x[1] = -x[1] ;
  x[2] = -x[2] ;
  x[3] = -x[3] ;
  #else
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) { x[mu] = -x[mu] ; }
  #endif
  //  contrary to what its name may suggest get_site_2piBZ translates
  //  the -pi to pi BZ coordinates to a site in the 0->2Pi lattice
  //  which is what the FFT uses
  return get_site_2piBZ( x , ND ) ;
}

#ifndef U1_DHT

// periodic fields using the DFT

static void 
periodic_dft( GLU_complex *__restrict *__restrict fields )
{
  // we know that there are no self-conjugate momenta beyond this point
  const size_t SYMM_POINT = compute_zeropoint( ) + 1 ;
  int *count = calloc( SYMM_POINT , sizeof( int ) ) ; // set up a counter

  size_t i ;
  // openmp does not play nice with the rng
#pragma omp parallel for private(i)
  for( i = 0 ; i < SYMM_POINT ; i++ ) {
    const uint32_t thread = get_GLU_thread( ) ;
    if( likely( count[i] == CONJUGATE_NOT_IN_LIST ) ) {
      count[i] = CONJUGATE_IN_LIST; // set the element of the list to 1
      const size_t b = conjugate_site( i ) ;
      size_t mu ;
      if( unlikely( i == b ) ) {
        #if ND%2 == 0
	for( mu = 0 ; mu < ND ; mu+=2 ) {
	  register const GLU_complex cache = par_polar_box( thread ) ;
	  fields[mu  ][i] = creal( cache ) ;
	  fields[mu+1][i] = cimag( cache ) ;
	}
        #else
	fields[0][i] = creal( par_polar_box( thread ) ) ;
	for( mu = 1 ; mu < ND ; mu+=2 ) {
	  register const GLU_complex cache = par_polar_box( thread ) ;
	  fields[mu  ][i] = creal( cache ) ;
	  fields[mu+1][i] = cimag( cache ) ;
	}
        #endif
      } else {
	for( mu = 0 ; mu < ND ; mu++ ) {
	  register const GLU_complex cache = par_polar_box( thread ) ;
	  fields[mu][i] = cache ;
	  fields[mu][b] = conj( cache ) ;
	}
      }
      // OK, so we can have conjugates that are less than the
      // symmetric 0 point and we have to accommodate for this
      if( b < SYMM_POINT ) { // # decreases with volume
	count[b] = CONJUGATE_IN_LIST ;
      }
      ////////////////////////////////////////
    }
  }
  free( count ) ;
  return ;
}

#else

// The DHT's periodicity requirements are basically the same, just use
// real fields ...
static void 
periodic_dht( GLU_real *__restrict *__restrict fields )
{
  // again we know that there are no self-conjugate momenta beyond this point
  //const int SYMM_POINT = LVOLUME ; was basically the old code
  const int SYMM_POINT = compute_zeropoint( ) + 1 ; 
  int *count = calloc( SYMM_POINT , sizeof(int) ) ; // set up a counter

  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < SYMM_POINT ; i++ ) {
    const uint32_t thread = get_GLU_thread( ) ;
    if( likely( count[i] == CONJUGATE_NOT_IN_LIST ) ) {
      count[i] = CONJUGATE_IN_LIST ; // set the element of the list to 1
      // get the momenta at "i" in the -pi to pi BZ 
      const size_t b = conjugate_site( i ) ;
      size_t mu ;
      #if ND%2 == 0
      for( mu = 0 ; mu < ND ; mu += 2 ) {
	register const GLU_complex cache = r2 * par_polar_box( thread ) ;
	fields[ mu ][ i ] = fields[ mu ][ b ] = creal( cache ) ;
	fields[ mu + 1 ][ i ] = fields[ mu + 1 ][ b ] = cimag( cache ) ;
      }
      #else
      fields[ 0 ][ i ] = fields[ 0 ][ b ] = r2 * creal( par_polar_box( thread ) ) ;
      for( mu = 1 ; mu < ND ; mu += 2 ) {
	register const GLU_complex cache = r2 * par_polar_box( thread ) ;
	fields[ mu ][ i ]     = fields[ mu ][ b ]     = creal( cache ) ;
	fields[ mu + 1 ][ i ] = fields[ mu + 1 ][ b ] = cimag( cache ) ;
      }
      #endif

      if( b < SYMM_POINT ) {
	count[b] = CONJUGATE_IN_LIST ; // set count[b] = 1 as well, it is already included
      }
      // and we are finished
    }
  }

  free( count ) ;
  return ;
}

#endif

// create the U1 fields
static int
create_u1( GLU_real *__restrict *__restrict U ,
	   const GLU_real alpha )
{
#ifndef HAVE_FFTW3_H

  fprintf( stderr , "[U(1)] Cannot create U(1) fields, "
	   "requires FFTW ... Leaving\n" ) ;
  return GLU_FAILURE ;

#else
  if( parallel_ffts( ) == GLU_FAILURE ) {
    fprintf( stderr , "[PAR] Problem with initialising the OpenMP "
	     "FFTW routines \n" ) ;
    return GLU_FAILURE ;
  }

  // alpha = 0.0795775387 is beta = 1 , gives noncompact plaq = 0.25 , test
  const GLU_real Nbeta = LVOLUME / ( ND * MPI * alpha ) ;
  fprintf( stdout , "\n[U(1)] 1/(%d Beta) :: %f \n\n" , ND , ( MPI * alpha ) ) ;

  size_t i ;

  // begin timing
  start_timer( ) ;

  // initialise the rng
  initialise_par_rng( NULL ) ;

#ifdef U1_DHT

  fprintf( stdout , "[U(1)] Using the DHT U1 code ... \n\n" ) ;

  GLU_real **in = fftw_malloc( ND * sizeof( GLU_real* ) ) ;
  GLU_real **out = fftw_malloc( ND * sizeof( GLU_real* ) ) ;
  for( i = 0 ; i < ND ; i++ ) {
    in[i] = ( GLU_real* )fftw_malloc( LVOLUME * sizeof( GLU_real ) ) ;
    out[i] = ( GLU_real* )fftw_malloc( LVOLUME * sizeof( GLU_real ) ) ;
  }

  fftw_plan plan[ ND ] ;
  create_plans_DHT( plan , in , out , Latt.dims , ND , ND ) ;

  periodic_dht( in ) ;

#else

  fprintf( stdout , "[U(1)] Using the DFT U1 code ... \n\n" ) ;

  GLU_complex **in = fftw_malloc( ND * sizeof( GLU_complex* ) ) ;
  GLU_complex **out = fftw_malloc( ND * sizeof( GLU_complex* ) ) ;
  for( i = 0 ; i < ND ; i++ ) {
    in[i] = ( GLU_complex* )fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;
    out[i] = ( GLU_complex* )fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;
  }
  
  fftw_plan plan[ ND ] , backward[ ND ] ;
  create_plans_DFT( plan , backward , in , out , Latt.dims , ND , ND ) ;

  periodic_dft( in ) ;

#endif

  // set up the distribution with the lattice mom
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    int flag = 0 ;
    const GLU_real f = 1.0 / (GLU_real)sqrt( gen_p_sq_feyn( i , &flag ) ) ;

    // these should become SIMD'd or something
    if( unlikely ( flag == 1 ) ) {
      #if ND == 4
      in[0][i] = in[1][i] = in[2][i] = in[3][i] = 0. + I * 0. ; 
      #else
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	in[mu][i] = 0. ;
      }
      #endif
    } else {
      #if ND == 4
      in[0][i] *= f ;
      in[1][i] *= f ;
      in[2][i] *= f ;
      in[3][i] *= f ;
      #else
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	in[mu][i] *= f ;
      }
      #endif
    }
  }

  // fft the fields allow for the parallel omp-ified fftws
#ifdef OMP_FFTW
  static size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    PSPAWN fftw_execute( plan[ mu ] ) ;
  }
  PSYNC ;
#else
  static size_t mu ;
#pragma omp parallel for private(mu) 
  for( mu = 0 ; mu < ND ; mu++ ) {
    PSPAWN fftw_execute( plan[ mu ] ) ;
  }
  PSYNC ;
#endif

  const GLU_real rbeta = 1.0 / (GLU_real)sqrt( Nbeta ) ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ )  {
    size_t nu ;
    for( nu = 0 ; nu < ND ; nu ++ ) {
      U[ nu ][ i ] = out[ nu ][ i ] * rbeta ;
    }
  }

  print_time( ) ;

  // FREE STUFF
 for( mu = 0 ; mu < ND ; mu++ ) {
   fftw_free( in[mu] ) ;
   fftw_free( out[mu] ) ;
   fftw_destroy_plan( plan[ mu ] ) ;  
   #ifndef U1_DHT // DHT's forward and backward plans are the same
   fftw_destroy_plan( backward[ mu ] ) ;  
   #endif
 }
 fftw_free( in ) ;
 fftw_free( out ) ;
 fftw_cleanup( ) ; 
#ifdef OMP_FFTW
 fftw_cleanup_threads( ) ; 
#endif

#endif
 return GLU_SUCCESS ;
}

#endif // HAVE_FFTW3_H

// Possibly a check that the plaquette has changed ..
int
suNC_cross_u1( struct site *__restrict lat , 
	       const struct u1_info U1INFO )
{
#ifndef HAVE_FFTW3_H
  fprintf( stderr , "[U1] Require FFTW to be linked to do quenched U(1)\n" ) ;
  return GLU_FAILURE ;
#else
  size_t i ; 
  GLU_real **U = malloc( ND * sizeof( GLU_real* ) ) ;
  for( i = 0 ; i < ND ; i++ ) {
    U[i] = ( GLU_real* )malloc( LVOLUME * sizeof( GLU_real ) ) ;
  }

  // create the U1 field ...
  create_u1( U , U1INFO.alpha ) ;
  
  // compute some U1 observables why not ?
  compute_U1_obs( (const GLU_real**)U , U1INFO.meas ) ;

  // multiply the < exponentiated > fields
#pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu ++ ) {
      // call to cexp is expensive
      register const GLU_real theta = U[mu][i] * U1INFO.charge ;
      const GLU_complex U1 = cos( theta ) + I * sin( theta ) ;
      size_t element ;
      for( element = 0 ; element < NCNC ; element++ ) {
	lat[i].O[mu][element] *= U1 ;
      }
    }
  }

  // free memory and stuff
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    free( U[mu] ) ;
  }
  free( U ) ;
#endif
  return GLU_SUCCESS ;
}
