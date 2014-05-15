/*
    Copyright 2013 Renwick James Hudspith

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

#include "geometry.h"     // general config-space geometry
#include "GLU_bswap.h"    // byte swaps
#include "cut_output.h"   // to automatically format our output file
#include "SM_wrap.h"      // for the smearing wrapper
#include "cut_routines.h" // for the mom cuts

// compute a short polyakov line
static void
small_poly( GLU_complex poly[ NCNC ] ,
	    const struct site *__restrict lat ,
	    const int site , const int dir , const int length )
{
  equiv( poly , lat[site].O[dir] ) ;
  // poly loop calculation ...
  int next = site , t ; 
  for( t = 1 ; t < length ; t++ ) {
    next = lat[next].neighbor[dir] ; 
    multab_atomic_right( poly , lat[next].O[dir] ) ;
  }
  return ;
}

// computes a polyakov line, looping around the dimension "dir"
double complex
poly( const struct site *__restrict lat , 
      int dir )
{
  // if you have stupidly set the dimension to something unreasonable
  // default the direction to ND
  if( dir > ND || dir < 0 ) { dir = ND ; }
  double complex sum = 0 ;
  int i ; 
#pragma omp parallel for private(i) reduction(+:sum)
  for( i = 0 ; i < LCU ; i++ ) {
    GLU_complex poly[ NCNC ] ;
    // use the correct site for one of the hypercubes ...
    int x[ ND ] ;
    get_mom_2piBZ( x , i , dir ) ;
    const int k = gen_site( x ) ;
    small_poly( poly , lat , k , dir , Latt.dims[dir] ) ;
    sum = sum + (double complex)trace( poly ) ;
  }
  return sum ;
}



// check if i and j form a triplet
#ifdef THREE_QUARK_CORRELATOR

static const double r3 = 1.7320508075688772 ;
static const double TWOPI_3 = TWOPI / 3. ;

static double
get_ang( const int V1[ ND-1 ] , const double V1_len , 
	 const int V2[ ND-1 ] , const double V2_len )
{
  if( V1_len == 0. || V2_len == 0. ) return M_PI ;
  int mu , scalar_prod = 0 ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    scalar_prod = V1[ mu ] * V2[ mu ] ;
  }
  return acos( scalar_prod / ( V1_len*V2_len ) ) ;
}

// compute LMINSQ
static double
compute_LMINSQ( const double a ,
		const double b ,
		const double c )
{
  return 0.5 * ( ( a*a + b*b + c*c ) +
		 r3 * sqrt( ( a+b+c ) * ( b-a+c ) *
			    ( a-b+c ) * ( a+b-c ) ) ) ;
}

static GLU_bool
is_a_triplet( const int i , 
	      const int j , 
	      double *rsq )
{
  // compute the vectors
  int VA[ ND-1 ] , VB[ ND-1 ] , VC[ ND-1 ] ; 
  get_mom_2piBZ( VA , i , ND-1 ) ;
  get_mom_2piBZ( VC , j , ND-1 ) ;
  int mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    VC[ mu ] = -VC[ mu ] ;
    VB[ mu ] = VA[ mu ] + VC[ mu ] ;
  }

  // compute their lengths
  double a = 0. , b = 0. , c = 0. ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    a += VA[ mu ] * VA[ mu ] ;
    b += VB[ mu ] * VB[ mu ] ;
    c += VC[ mu ] * VC[ mu ] ;
  }
  a = sqrt( a ) ; b = sqrt( b ) ; c = sqrt( c ) ; 

  // compute the angle between these ...
  if( get_ang( VA , a , VB , b ) < TWOPI_3 ) return GLU_FALSE ;
  else if( get_ang( VA , a , VC , c ) < TWOPI_3 ) return GLU_FALSE ;
  else if( get_ang( VB , b , VC , c ) < TWOPI_3 ) return GLU_FALSE ;

  // compute the LSQ and return TRUE
  *rsq = compute_LMINSQ( a , b , c ) ;

  return GLU_TRUE ;
}
#endif

/// correlator
static void
static_quark_correlator( double *__restrict result ,
			 double *__restrict trtr ,
			 const GLU_complex *__restrict *__restrict poly ,
			 const struct veclist *__restrict list ,
			 const int rsq_count )

{
  int i ;
  // now we can loop the triplet computing the correlators ...
#pragma omp parallel for private(i)
  for( i = 0 ; i < rsq_count ; i++ ) {

    // COOL, can now contract these over the volume
    GLU_complex res ;
    register double tr = 0.0 ;
    register double tracetrace = 0.0 ;
    const int ref_idx = list[i].idx ;

    int t ;
    for( t = 0 ; t < Latt.dims[ ND-1 ] ; t++ ) {
      int k ;
      for( k = 0 ; k < LCU ; k++ ) {
	const int idx1 = ( k + LCU * t ) ;
	const int idx2 = ( k + ref_idx ) < LCU ? ( ref_idx + idx1 ) : ( ref_idx + idx1 - LCU ) ;

	// really need to think hard about this one
	#ifdef THREE_QUARK_CORRELATOR
	const int idx3 = ( i + list[k].idx_2 )%LCU + LCU * t ;

	// extract the singlet
	double complex trtrtr = (double complex)( trace( poly[idx1] ) * 
						  trace( poly[idx2] ) * 
						  trace( poly[idx3] ) ) ;

	// double trace terms
	trace_ab( &res , poly[idx2] , poly[idx3] ) ;
	double complex TrPxTrPyPz = -9.0 * trace( poly[idx1] ) * res ;
	trace_ab( &res , poly[idx1] , poly[idx3] ) ;
	double complex TrPyTrPxPz = -9.0 * trace( poly[idx2] ) * res ;
	trace_ab( &res , poly[idx1] , poly[idx2] ) ;
	double complex TrPzTrPxPy = -9.0 * trace( poly[idx3] ) * res ;

	// and the triple trace terms
	double complex TrPxPyPz ;
	trace_abc( &TrPxPyPz , poly[idx1] , poly[idx2] , poly[idx3] ) ;
	double complex TrPxPzPy ;
	trace_abc( &TrPxPzPy , poly[idx1] , poly[idx3] , poly[idx2] ) ;

	// this is the singlet
	tr += creal( +27. * ( trtrtr )
		     -9.0 * ( TrPxTrPyPz + TrPyTrPxPz + TrPzTrPxPy )
		     +3.0 * ( TrPxPyPz + TrPxPzPy ) ) / 6.0 ;

	// the decuplet
	tracetrace += creal( +27. * trtrtr 
			     +9.0 * ( TrPxTrPyPz + TrPyTrPxPz + TrPzTrPxPy )
			     +3.0 * ( TrPxPyPz + TrPxPzPy ) ) ) / 6.0 ;
	  #if 0
	  //Here I list the other possible terms
	  
	  // the decuplet
	  result[i] += ( +27. * trtrtr 
	  +9.0 * ( TrPxTrPyPz + TrPyTrPxPz + TrPzTrPxPy )
	  +3.0 * ( TrPxPyPz + TrPxPzPy ) ) / 6.0 ;

	  // the octet
	  result[i] += ( +27. * trtrtr 
	  +9.0 * ( TrPxTrPyPz - TrPzTrPxPy )
	  +3.0 * ( TrPxPyPz ) ) / 6.0 ;

	  // the octet'
	  result[i] += ( +27. * trtrtr 
	  -9.0 * ( TrPxTrPyPz - TrPzTrPxPy )
	  +3.0 * ( TrPxPzPy ) ) / 6.0 ;
	  #endif

	#else

	// trace of the product
	trace_ab_dag( &res , poly[idx1] , poly[idx2] ) ;
	tr += (double)creal( res ) ;

	// and the trace-trace
	const GLU_complex tr1 = trace( poly[idx1] ) ;
	const GLU_complex tr2 = trace( poly[idx2] ) ;
	tracetrace += creal( tr1 ) * creal( tr2 ) + cimag( tr1 ) * cimag( tr2 ) ;

	#endif
      }
     }
    result[i] = tr ;
    trtr[i]   = tracetrace ;
  }
  return ;
}

// and so for the correlator measurements
void
Coul_staticpot( struct site *__restrict lat ,
		const struct cut_info CUTINFO ,
		const struct sm_info SMINFO )
{
  // do some smearing? Overwrites lat
  SM_wrap_struct( lat , SMINFO ) ;

  // important!! T is the length of the polyakov loop in time direction
  const int T = CUTINFO.max_t ;

  printf( "\n[STATIC-POTENTIAL] measurements at T = %d\n" , T ) ;
#ifndef THREE_QUARK_CORRELATOR
  printf( "[STATIC-POTENTIAL] cut off at r^2 = %d\n" , CUTINFO.max_mom ) ;
#endif

  // compute the ratios of the dimensions in terms of the smallest
  simorb_ratios( ND ) ;

  // compute the momentum list
  int size[1] = {} ;
  struct veclist *list = compute_veclist( size , CUTINFO , ND-1 , GLU_TRUE ) ;

  printf( "SIZE :: %d \n" , size[0] ) ;

  // precompute the correlator
  // compute all of the poly loops, lattice-wide
  GLU_complex **poly = malloc( LVOLUME * sizeof( GLU_complex* ) ) ;
  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    poly[i] = (GLU_complex*) malloc( NCNC * sizeof( GLU_complex ) ) ;
    identity( poly[i] ) ;
  }

  // set up the outputs
  const char *str = output_str_struct( CUTINFO ) ;  
  FILE *Ap = fopen( str , "wb" ) ; 

  // write the list, in cut_outputs.c
  write_mom_veclist( Ap , size , list , ND-1 ) ;

  // tell us where it is going
  printf( "[STATIC-POTENTIAL] writing correlators up to T = %d separation \n" , T ) ;
  printf( "[STATIC-POTENTIAL] writing output file to %s \n" , str ) ;

  // write out the size of the timeslices we are looking at
  int Tcorrs[1] = { T-1 } ;
  if( !WORDS_BIGENDIAN ) { bswap_32( 1 , Tcorrs ) ; }
  fwrite( Tcorrs , sizeof(int) , 1 , Ap ) ;
  if( !WORDS_BIGENDIAN ) { bswap_32( 1 , Tcorrs ) ; }  

  // allocate the results
  double *result = malloc( size[0] * sizeof( double ) ) ; 
  double *trtr   = malloc( size[0] * sizeof( double ) ) ; 

  // precompute the polyakov loop
  int t ; // t == 0 term makes no sense ...
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
      GLU_complex temp[ NCNC ] ;
      equiv( temp , poly[i] ) ;
      const int tpos = ( i + t * LCU ) % LVOLUME ;
      multab_suNC( poly[i] , temp , lat[tpos].O[ND-1] ) ;
    }
    
    // compute the three quark correlator
    static_quark_correlator( result , trtr , (const GLU_complex**)poly , 
			     list , size[0] ) ;

    // write out the result of all that work
    write_g2g3_to_list( Ap , result , trtr , size ) ;

    printf( "[STATIC-POTENTIAL] T = %d sepation computed and written \n" , t ) ;
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

  return ;
}

