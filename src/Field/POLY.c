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
  PFOR( i = 0 ; i < rsq_count ; i++ ) {

    // COOL, can now contract these over the volume
    GLU_real res ;
    register double tr = 0.0 ;
    register double tracetrace = 0.0 ;

    // the (positive) lattice vector for the separation
    int separation[ ND ] ;
    get_mom_2piBZ( separation , list[i].idx , ND-1 ) ;

    // loop the lattice varying source and sink with the correct separation
    int k = 0 ;
    for( k = 0 ; k < LCU ; k++ ) {
      
      // translate the source from k by a vector separation
      const int translate = compute_spacing( separation , k , ND - 1 ) ;
      
      int t = 0 ;
      for( t = 0 ; t < Latt.dims[ ND-1 ] ; t++ ) {
	
	// compute the source and sink for this time separation
	const int source = k + LCU * t ;
	const int sink   = translate + LCU * t ;
	
	// trace of the product
	trace_ab_dag_Re( &res , poly[source] , poly[sink] ) ;
	tr += (double)res ;
	
	// and the trace-trace
	const GLU_complex tr1 = trace( poly[source] ) ;
	const GLU_complex tr2 = trace( poly[sink] ) ;
	tracetrace += creal( tr1 ) * creal( tr2 ) + cimag( tr1 ) * cimag( tr2 ) ;
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
  // do some smearing? Overwrites lat, should I state that I only think
  // SPATIAL dimensional smearing is safe ... probably
  SM_wrap_struct( lat , SMINFO ) ;

  // important!! T is the length of the polyakov loop in time direction
  const int T = CUTINFO.max_t ;

  //
  printf( "\n[STATIC-POTENTIAL] measurements at T = %d\n" , T ) ;

  // compute the ratios of the dimensions in terms of the smallest
  simorb_ratios( ND ) ;

  // compute the momentum list
  int size[1] = {} ;
  struct veclist *list = compute_veclist( size , CUTINFO , ND-1 , GLU_TRUE ) ;

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
      multab_atomic_right( poly[i] , lat[( i + t * LCU ) % LVOLUME].O[ND-1] ) ;
    }
    
    // compute the three quark correlator?
    static_quark_correlator( result , trtr , (const GLU_complex**)poly , 
			     list , size[0] ) ;

    // write out the result of all that work
    write_g2g3_to_list( Ap , result , trtr , size ) ;

    // and tell us which temporal separation has been done
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

