/*
Copyright 2013-2025 Renwick James Hudspith

    This file (MOMggg.c) is part of GLU.

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
   @file MOMggg.c
   @brief computation of the non-exceptional three point function and 
   the gluon propagator
 */
#include "Mainfile.h"

#include "cut_output.h"  // output file formatting
#include "geometry.h"    // general geometry look-ups
#include "lie_mats.h"    // lie algebra stuff, fs and ds
#include "triplet_gen.h" // equivalent triplet calculator

//#define LIE_PROJECTION

#ifdef VIEW_TRIPLETS
static void
print_momenta( int p1[ND] , int p2[ND] )
{
  size_t mu ;
  fprintf( stdout , "(" ) ; for( mu = 0 ; mu < ND ; mu++ ) { 
    fprintf( stdout , "%d " , p1[mu] ) ; } ; fprintf( stdout , ") " ) ;
  printf( "(" ) ; for( mu = 0 ; mu < ND ; mu++ ) { 
    fprintf( stdout , "%d " , p2[mu] ) ; } ; fprintf( stdout , ") " ) ;
  fprintf( "(" ) ; for( mu = 0 ; mu < ND ; mu++ ) { 
    fprintf( stdout , "%d " , -(p1[mu]+p2[mu]) ) ; } ; 
  fprintf( stdout , ") \n" ) ;
}
#endif

//// COMPUTES THE NONEXCEPTIONAL MOMENTUM CONFIGURATION 3POINT FUNCTION
int
write_nonexceptional_g2g3( FILE *__restrict Ap , 
			   const struct site *__restrict A , 
			   const struct veclist *__restrict list , 
			   size_t num_mom[ 1 ] , 
			   const size_t nnmax )
{
#ifdef LIE_PROJECTION
  init_generators() ;
  compute_fs_and_ds() ;
#endif

  const double g2_norm = 2.0 / ( 3.0 * ( NCNC - 1 ) * ( ND - 1 ) * LVOLUME ) ;
  // factor of 3 comes from the average over triplets
#ifdef LIE_PROJECTION
  const double g3_norm = 1.0 / ( LVOLUME * NC * ( NCNC - 1 ) ) ;
#else
  const double g3_norm = 4.0 / ( LVOLUME * NC * ( NCNC - 1 ) ) ;
#endif
 
  //need to compute the triplet, first we compute the momentum of our list
  int **momentum = malloc( num_mom[ 0 ] * sizeof( int * ) ) ; 

  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < num_mom[0] ; i++ ) {
    momentum[ i ] = ( int * )malloc( ND * sizeof ( int ) ) ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      momentum[ i ][ mu ] = list[i].MOM[mu] ;
    }
  }

  // fill up trip by calling get_triplet again .. read from a file ..
  // set up the number of triplets ... 
  int *trip = malloc( nnmax/2 * sizeof( int ) ) ;

  // allocate both the triplet and the projector
  int **triplet = NULL ;
  double **proj = NULL ;

  // malloc these as they are pretty big ...
  double *g2 = NULL , *g3 = NULL ;

  int flag = GLU_SUCCESS ;
  size_t checker = 0 ;

  // look for a file or just calculate it the dumb way
  read_trip( trip , nnmax ) ;

  size_t nn ;
  size_t counter[ 1 ] = { 1 } ;
  for( nn = 0 ; nn < nnmax/2 ; nn++ ) {
    counter[ 0 ] += trip[ nn ] ;
  }
  
  // set this up so that it will never be altered, important.
  const size_t count = counter[ 0 ]  ;
  if( count == 0 ) { 
    flag = GLU_FAILURE ;
    goto memfree ;
  }

  triplet = malloc( count * sizeof( int* ) ) ;
  proj = malloc( count * sizeof( double* ) ) ;

#pragma omp parallel for private(i)
  for( i = 0 ; i < count ; i++ ) {
    proj[i] = ( double* )malloc( ND * ND * ND *  sizeof ( double ) ) ;
    triplet[i] = ( int* )malloc( 3 *  sizeof( int ) ) ;
  }

  g2 = malloc( count * sizeof( double ) ) ;
  g3 = malloc( count * sizeof( double ) ) ;

  // read in the triplet and the projector 
  flag = read_triplet_and_proj( triplet , // reads in the trips
				proj ,  // computes the projector
				(const int**)momentum , // momlist
				nnmax , // maximum momentum
				count , //size of the trip array
				num_mom[ 0 ] ) ; // size mom list
  
  if( flag == GLU_FAILURE ) {
    goto memfree ;
  }
  
  write_triplet_mom_list( Ap , counter , momentum , triplet ) ;

  for( nn = 0 ; nn < nnmax/2 ; nn ++ ) {

    double psq = 0. ;
    double oneOpsq = 0. ;
    if( trip[ nn ] != 0 ) {
      // compute psq
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	const double cache = momentum[ triplet[ checker ][ 0 ] ][ mu ] * \
	  Latt.twiddles[ mu ] ; 
        #ifdef PSQ_MOM
	const double temp =  cache ;
        #else 
	// Sin variant hmmm, psq's may (will) be different in this case
	const double temp =  2. * sin( 0.5 * cache );
        #endif
	psq += temp * temp ;
      }
	
      if( psq < PREC_TOL ) {
	psq = 1.0 ;
      }

      // wrap the psq into the g3 normalisation ...
      oneOpsq = g3_norm / ( psq ) ;
    }

    // count through the available triplets .. should be able to parallelize
    for( i = 0 ; i < (size_t)trip[ nn ] ; i++ ) { 	  
      const size_t place = checker + i ;
	  
      // initialise the threepoint
      g3[ place ] = 0. ;
      
      const size_t trip0 = list[ triplet[ place ][ 0 ] ].idx ; 
      const size_t trip1 = list[ triplet[ place ][ 1 ] ].idx ; 
      const size_t trip2 = list[ triplet[ place ][ 2 ] ].idx ; 

      // apply projector ? try this one 
      size_t z ;
      double sum = 0. ;
      #pragma omp parallel for private(z) reduction(+:sum)
      for( z = 0 ; z < ( ND * ND * ND ) ; z++ ) {
	if( fabs( proj[ place ][ z ] ) > PREC_TOL ) {
	  // faster , simpler triplet multiplication
	  GLU_complex test = 0. ;
	  const size_t rho = z % ND ;
	  const size_t nu = ( ( z - z % ND ) / ND ) % ND ;
	  const size_t mu = ( ( z - z % ( ND * ND ) ) / ( ND * ND ) ) % ND ;
	  
	  // do the projection either in lie space or not
          #ifdef LIE_PROJECTION
	  test = (GLU_complex)ifabc_ABC( A[ trip0 ].O[mu] , 
					 A[ trip1 ].O[nu] , 
					 A[ trip2 ].O[rho] ) ;
          #else
	  trace_abc( &test ,
		     A[ trip0 ].O[mu] ,
		     A[ trip1 ].O[nu] ,
		     A[ trip2 ].O[rho] ) ;
          #endif
	  
	  // and sum it ... Imaginary part cancels !
          #ifdef CUT_FORWARD
	  sum = sum + (double)( proj[ place ][ z ] * (double)creal( test ) ) ;
          #else
	  sum = sum - (double)( proj[ place ][ z ] * (double)creal( test ) ) ;
          #endif
	}
      }
      
      // psq factor implicit in the norm, g3_norm included too
      g3[ place ] = sum * oneOpsq ;
      
      // compute the gluonic two point function here ...
      double res = 0. ;
      size_t mu ;
      #pragma omp parallel for private(mu) reduction(+:res)
      for( mu = 0 ; mu < ND ; mu++ ) {
	GLU_complex tr ;
	trace_ab_dag( &tr , A[ trip0 ].O[mu] , A[ trip0 ].O[mu] ) ;
	register double loc_sum = (double)creal( tr ) ;
	trace_ab_dag( &tr , A[ trip1 ].O[mu] , A[ trip1 ].O[mu] ) ;
	loc_sum += (double)creal( tr ) ;
	trace_ab_dag( &tr , A[ trip2 ].O[mu] , A[ trip2 ].O[mu] ) ;
	loc_sum += (double)creal( tr ) ;	
	res = res + (double)loc_sum ;
      }
      g2[ place ] = res * g2_norm ; //three contributing triplets
    }
    // a value for which we add "place" to in the files..
    checker += trip[nn] ;
  }
  
  // we write out the two and three point functions in the output file
  write_g2g3_to_list( Ap , g2 , g3 , counter ) ; 

  // great success!
  flag = GLU_SUCCESS ;

 memfree :

  // free the array denoting the size of triplets
  free( trip ) ;

  // free the two and three point functions
  free( g2 ) ;
  free( g3 ) ;

#pragma omp parallel for private(i)
  for( i = 0 ; i < num_mom[ 0 ] ; i++ ) {
    free( momentum[ i ] ) ;
  }
  free( momentum ) ;

#pragma omp parallel for private(i)
  for( i = 0 ; i < count ; i++ ) {
    if( triplet != NULL ) {
      free( triplet[ i ] ) ;
    }
    if( proj != NULL ) {
      free( proj[i] ) ;
    }
  }
  free( triplet ) ;
  free( proj ) ;

#ifdef LIE_PROJECTION
  free_f_and_d( ) ;
  free_generators( ) ;
#endif

  return flag ;
}

#ifdef LIE_PROJECTION
 #undef LIE_PROJECTION
#endif
