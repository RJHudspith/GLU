/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (POLY.c) is part of GLU.

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
   @file POLY.c
   @brief Coulomb gauge fixed polyakov loop correlators 

   I write out the r^2's and the real part of the singlet and q-q correlators

   @warning Calls the smearing wrapper, hence overwrites the lattice links
 */

#include "Mainfile.h"

#include "cut_output.h"   // to automatically format our output file
#include "cut_routines.h" // for the mom cuts
#include "geometry.h"     // general config-space geometry
#include "GLU_bswap.h"    // byte swaps
#include "plan_ffts.h"    // convolution
#include "SM_wrap.h"      // for the smearing wrapper

// compute a short polyakov line
static void
small_poly( GLU_complex poly[ NCNC ] ,
	    const struct site *__restrict lat ,
	    const size_t site , 
	    const size_t dir , 
	    const size_t length )
{
  equiv( poly , lat[site].O[dir] ) ;
  // poly loop calculation ...
  size_t next = site , t ; 
  for( t = 1 ; t < length ; t++ ) {
    next = lat[next].neighbor[dir] ; 
    multab_atomic_right( poly , lat[next].O[dir] ) ;
  }
  return ;
}

// computes a polyakov line, looping around the dimension "dir"
double complex
poly( const struct site *__restrict lat , 
      size_t dir )
{
  // if you have stupidly set the dimension to something unreasonable
  // default the direction to ND
  if( dir > ND ) { dir = ND ; }
  double complex sum = 0 ;
  size_t i , subvolume = 1 ;
  for( i = 0 ; i < ND ; i++ ) {
    subvolume *= (i!=dir)? Latt.dims[i] : 1 ;
  }
#pragma omp parallel for private(i) reduction(+:sum)
  for( i = 0 ; i < subvolume ; i++ ) {
    GLU_complex poly[ NCNC ] ;
    // use the correct site for one of the hypercubes ...
    int x[ ND ] ;
    get_mom_2piBZ( x , i , dir ) ;
    const size_t k = gen_site( x ) ;
    small_poly( poly , lat , k , dir , Latt.dims[dir] ) ;
    sum = sum + (double complex)trace( poly ) ;
  }
  return sum ;
}

// If we have FFTW we use it for the convolutions instead of our slow
// method, which is slow but equivalent
#ifdef HAVE_FFTW3_H

// trace-trace computation is nice and easy
static void
compute_trtr( double complex *__restrict trtr ,      
	      GLU_complex *__restrict out ,
	      GLU_complex *__restrict in ,
	      const GLU_complex *__restrict *__restrict poly ,
	      const struct veclist *__restrict list ,
	      const fftw_plan forward ,
	      const fftw_plan backward ,
	      const size_t rsq_count ,
	      const size_t t )
{
  size_t i ;
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i++ ) {
    in[ i ] = trace( poly[ i + LCU*t ] ) ;
  }
  // FFT
  fftw_execute( forward ) ;
  // convolve
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i++ ) {
    out[i] *= conj( out[i] ) ;
  }
  // FFT
  fftw_execute( backward ) ;
  // and set the trace-trace
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < rsq_count ; i++ ) {
    trtr[ i ] += in[ list[i].idx ] ;
  }
  return ;
}

// traced computation is hard, I want to do the computation in momentum space
// but cannot overwrite poly?
static void
compute_tr( double complex *__restrict tr ,      
	    GLU_complex *__restrict out ,
	    GLU_complex *__restrict in ,
	    GLU_complex *__restrict *__restrict slice_poly ,
	    const GLU_complex *__restrict *__restrict poly ,
	    const struct veclist *__restrict list ,
	    const fftw_plan forward ,
	    const fftw_plan backward ,
	    const size_t rsq_count ,
	    const size_t t )
{
  // index by index FFTs into a temporary
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    equiv( slice_poly[ i ] , poly[ i + LCU * t ] ) ; 
  }
  // FFT all indices
  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    // set in 
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i++ ) {
      in[ i ] = slice_poly[ i ][ j ] ;
    }
    fftw_execute( forward ) ;
    // store the result in slice_poly again
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i++ ) {
      slice_poly[ i ][ j ] = out[ i ] ;
    }
  }
  // now we can contract
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i++ ) {
    trace_ab_dag( &out[ i ] , slice_poly[ i ] , slice_poly[ i ] ) ;
  }
  fftw_execute( backward ) ;
  // and set tr
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < rsq_count ; i++ ) {
    tr[ i ] += in[ list[ i ].idx ] ;
  }
  return ;
}

/// correlator
static int
static_quark_correlator( double complex *__restrict result ,
			 double complex *__restrict trtr ,
			 const GLU_complex *__restrict *__restrict poly ,
			 const struct veclist *__restrict list ,
			 const size_t rsq_count )

{
  // init parallel threads, maybe
  if( parallel_ffts( ) == GLU_FAILURE ) {
    fprintf( stderr , "[PAR] Problem with initialising the OpenMP "
	              "FFTW routines \n" ) ;
    return GLU_FAILURE ;
  }

  GLU_complex **slice_poly = NULL ;
  if( GLU_malloc( (void**)&slice_poly , 16 , LCU * sizeof( GLU_complex* ) ) 
      != 0 ) {
    return GLU_FAILURE ;
  }
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    GLU_malloc( (void**)&slice_poly[i] , 16 , NCNC * sizeof( GLU_complex ) ) ;
  }

  // FFTW routines
  GLU_complex *in = fftw_malloc( LCU * sizeof( GLU_complex ) ) ;
  GLU_complex *out = fftw_malloc( LCU * sizeof( GLU_complex ) ) ;

  // create some plans
  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , in , out , ND-1 ) ;

  size_t t ;
  for( t = 0 ; t < Latt.dims[ ND - 1 ] ; t++ ) {
    // compute the tr tr combination
    compute_trtr( trtr , out , in , poly , list , 
		  forward , backward , rsq_count , t ) ;

    // compute the traced
    compute_tr( result , out , in , slice_poly , poly , list , 
		forward , backward , rsq_count , t ) ;
  }

  // and normalise, LCU is the norm from the FFTs
  const double NORM = 1.0 / ( (double)LCU * Latt.dims[ ND-1 ] ) ;

  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < rsq_count ; i++ ) {
    result[ i ] *= NORM ;
    trtr[ i ] *= NORM ;
  }

  // free the FFTs
  fftw_destroy_plan( backward ) ;
  fftw_destroy_plan( forward ) ;
  fftw_cleanup( ) ;
#ifdef OMP_FFTW
  fftw_cleanup_threads( ) ;
#endif
  fftw_free( in ) ;
  fftw_free( out ) ;

  // free the temporary polyakov loops
#pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    free( slice_poly[ i ] ) ;
  }
  free( slice_poly ) ;

  return GLU_SUCCESS ;
}

// NON-fftw version is pretty slow
#else

/// correlator
static int
static_quark_correlator( double complex *__restrict result ,
			 double complex *__restrict trtr ,
			 const GLU_complex *__restrict *__restrict poly ,
			 const struct veclist *__restrict list ,
			 const size_t rsq_count )

{
  const double NORM = 1.0 / (double)Latt.dims[ ND - 1 ] ;
  size_t i ;
  // now we can loop the triplet computing the correlators ...
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < rsq_count ; i++ ) {

    // COOL, can now contract these over the volume
    GLU_complex res ;
    register double complex tr = 0.0 ;
    register double complex tracetrace = 0.0 ;

    // the (positive) lattice vector for the separation
    int separation[ ND ] ;
    get_mom_2piBZ( separation , list[i].idx , ND-1 ) ;

    // loop the lattice varying source and sink with the correct separation
    size_t k = 0 ;
    for( k = 0 ; k < LCU ; k++ ) {
      
      // translate the source from k by a vector separation
      const int translate = compute_spacing( separation , k , ND - 1 ) ;
      
      size_t t = 0 ;
      for( t = 0 ; t < Latt.dims[ ND-1 ] ; t++ ) {
	
	// compute the source and sink for this time separation
	const size_t source = k + LCU * t ;
	const size_t sink   = translate + LCU * t ;
	
	// trace of the product
	trace_ab_dag( &res , poly[source] , poly[sink] ) ;
	tr += res ;
	
	// and the trace-trace
	const GLU_complex tr1 = trace( poly[source] ) ;
	const GLU_complex tr2 = trace( poly[sink] ) ;
	tracetrace += tr1 * tr2 ;
      }
    }
    result[i] = tr * NORM ;
    trtr[i]   = tracetrace * NORM ;
  }

  return GLU_SUCCESS ;
}
#endif

// and so for the correlator measurements
int
Coul_staticpot( struct site *__restrict lat ,
		const struct cut_info CUTINFO ,
		const struct sm_info SMINFO )
{
  // do some smearing? Overwrites lat, should I state that I only think
  // SPATIAL dimensional smearing is safe ... probably
  SM_wrap_struct( lat , SMINFO ) ;

  // important!! T is the length of the polyakov loop in time direction
  const size_t T = CUTINFO.max_t ;

  // leave if the loop arguments make no sense
  if( T < 1 || LVOLUME < 1 ) return GLU_FAILURE ;

  fprintf( stdout , "\n[STATIC-POTENTIAL] measurements at T = %zu\n" , T ) ;

  // compute the ratios of the dimensions in terms of the smallest
  simorb_ratios( ND ) ;

  // precompute the correlator
  // compute all of the poly loops, lattice-wide
  GLU_complex **poly = NULL ;
  if( GLU_malloc( (void**)&poly , 16 , LVOLUME * sizeof( GLU_complex* ) ) != 0 ) {
    return GLU_FAILURE ;
  }

  size_t i ;
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    GLU_malloc( (void**)&poly[i] , 16 , NCNC * sizeof( GLU_complex ) ) ;
    identity( poly[i] ) ;
  }

  // compute the momentum list
  int size[1] = { 0 } ;
  struct veclist *list = compute_veclist( size , CUTINFO , ND-1 , GLU_TRUE ) ;

  // set up the outputs
  const char *str = output_str_struct( CUTINFO ) ;  
  FILE *Ap = fopen( str , "wb" ) ; 

  // write the list, in cut_outputs.c
  write_mom_veclist( Ap , size , list , ND-1 ) ;

  // tell us where it is going
  fprintf( stdout , "[STATIC-POTENTIAL] writing correlators up to "
	   "T = %zu separation \n" , T ) ;
  fprintf( stdout , "[STATIC-POTENTIAL] writing output file to %s \n" , str ) ;

  // write out the size of the timeslices we are looking at
  size_t Tcorrs[1] = { T-1 } ;
  if( !WORDS_BIGENDIAN ) { bswap_32( 1 , Tcorrs ) ; }
  fwrite( Tcorrs , sizeof(int) , 1 , Ap ) ;
  if( !WORDS_BIGENDIAN ) { bswap_32( 1 , Tcorrs ) ; }  

  // allocate the results
  double complex *result = malloc( size[0] * sizeof( double complex ) ) ; 
  double complex *trtr   = malloc( size[0] * sizeof( double complex ) ) ; 

  // precompute the polyakov loop
  size_t t ;
  for( t = 1 ; t < T ; t++ ) {

    // initialise results 
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < size[0] ; i++ ) {
      result[i] = trtr[i] = 0.0 ;
    }

    // atomically multiply link matrix at posit (i,t)
    // into the array of matrices poly
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
      multab_atomic_right( poly[i] , lat[( i + t * LCU ) % LVOLUME].O[ND-1] ) ;
    }
    
    // compute the two quark correlator
    static_quark_correlator( result , trtr , 
			     (const GLU_complex**)poly , 
			     list , size[0] ) ;

    // write out the result of all that work
    write_complex_g2g3_to_list( Ap , result , trtr , size ) ;

    // and tell us which temporal separation has been done
    fprintf( stdout , "[STATIC-POTENTIAL] T = %zu sepation computed "
	     "and written \n" , t ) ;
  }

  // close the file
  fclose( Ap ) ; 

  // memory free ...
  free( result ) ;
  free( trtr ) ;
  free( (void*)str ) ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) { 
    free( poly[i] ) ;
  }
  free( poly ) ;

  // free the momentum list
  free( list ) ;

  return GLU_SUCCESS ;
}
