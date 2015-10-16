/*
    Copyright 2013 Renwick James Hudspith

    This file (OrLandau.c) is part of GLU.

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
   @file OrLandau.c
   @brief Over relaxed Landau and Coulomb gauge fixing codes
 */

#include "Mainfile.h"       // general includes

#include "geometry.h"       // geom for the draughtboarding
#include "givens.h"         // su(2) rotations
#include "gftests.h"        // theta test
#include "gtrans.h"
#include "plaqs_links.h"    // plaquettes
#include "random_config.h"  // latt reunitisation

// inits for the draughtboard
static int *redsites , *blacksites ;
static int NRED = 0 , NBLACK = 0 ;

// free the draughtboard
static void
free_cb( void ) 
{
  free( redsites ) ;
  free( blacksites ) ;
}

// initialise the draughtboarding
static void
init_cb( const int LENGTH ) 
{
  int i , n[ ND ] ;
  for( i = 0 ; i < LENGTH ; i++ ) {
    get_mom_2piBZ( n , i , ND ) ;
    int mu , mode_sum = 0 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      mode_sum += n[ mu ] ;
    }
    if( mode_sum%2 == 0 ) {
      NRED++ ;
    } else {
      NBLACK++ ;
    }
  }
  // malloc and set
  redsites = malloc( NRED * sizeof( int ) ) ;
  blacksites = malloc( NBLACK * sizeof( int ) ) ;
  // set back to zero
  NRED = NBLACK = 0 ;
  for( i = 0 ; i < LENGTH ; i++ ) {
    get_mom_2piBZ( n , i , ND ) ;
    int mu , mode_sum = 0 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      mode_sum += n[ mu ] ;
    }
    if( mode_sum%2 == 0 ) {
      redsites[NRED] = i ;
      NRED++ ;
    } else {
      blacksites[NBLACK] = i ;
      NBLACK ++ ;
    }
  }
  return ;
}

// a single iteration of the overrelaxed gauge fixing
static void
OR_single( struct site *__restrict lat ,
	   const int su2_index ,
	   const double OrParam ,
	   const int i ,
	   const int DIMS )
{
  GLU_complex L[ NCNC ] ;
  zero_mat( L ) ;
  int mu , j , k ;
  // loop directions summing into L
  for( mu = 0 ; mu < DIMS ; mu++ ) {
    // compute U(x+\mu/2) + U^{dagger}(x-\mu/2)
    const int back = lat[i].back[mu] ;
    for( j = 0 ; j < NC ; j++ ) {
      for( k = 0 ; k < NC ; k++ ) {
	L[ k + j * NC ] += conj( lat[i].O[mu][j+k*NC] ) +
	  lat[back].O[mu][k+j*NC] ;
      }
    }
  }
  // hits the link to the left and the one to the right with
  // gauge transformation matrices
  GLU_complex s0 , s1 ;
  OrRotation( L , &s0 , &s1 , OrParam , su2_index ) ;

  // gauge rotate
  for( mu = 0 ; mu < ND ; mu++ ) {
    const int back = lat[i].back[mu] ;
    shortened_su2_multiply( lat[i].O[mu] , s0 , s1 , 
			    -conj(s1) , conj(s0) , su2_index ) ;
    shortened_su2_multiply_dag( lat[back].O[mu] , s0 , s1 , 
				-conj(s1) , conj(s0) , su2_index ) ;
  }

  return ;
}

// perform one iteration of the overrelaxed gauge fixing routine
static void
OR_iteration( struct site *__restrict lat ,
	      const int su2_index ,
	      const double OrParam ,
	      const int t ,
	      const int DIMS )
{
  // perform an overrelaxation step
  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < NRED ; i++ ) {  
    OR_single( lat , su2_index  , OrParam , 
	       redsites[i] + LCU * t , DIMS ) ;
  }
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < NBLACK ; i++ ) {    
    OR_single( lat , su2_index , OrParam , 
	       blacksites[i] + LCU * t , DIMS ) ;
  }
  return ;
}

// output the data, pass lat for the plaquette
static void
output_fixing_info( lat , theta , iters )
     struct site *__restrict lat ;
     const double theta ;
     const int iters ;
{
  // reunitarise just to limit the damage from round-off
  latt_reunitU( lat ) ;

  ////////// Print out the Gauge Fixing information /////////////

  printf( "[GF] Plaquette :: %1.15f \n[GF] Accuracy :: %1.4e\n" , 
	  av_plaquette( lat ) , theta ) ;
  GLU_real tr ;
  const double link = indivlinks( lat , &tr ) ;
  printf( "[GF] Iters :: %d\n[GF] Link trace :: %1.15f || Maximum :: %1.15f\n" , iters , link , tr / NC ) ; 
  double lin , log ;
  const_time( lat , &lin , &log ) ; 
  printf( "[GF] Temporal constance || Lin %e || Log %e \n" , lin , log ) ;
  gauge_functional( lat ) ;
  printf( "[GF] Functional :: %1.15f\n" , gauge_functional( lat ) ) ;

  ///////////////////////////////////////////////////////////////
  return ;
}

// gauge fix
int
OrLandau( struct site *__restrict lat ,
	  double *theta ,
	  const int MAX_ITERS , 
	  const double ACC ,
	  const double OrParam )
{
  GLU_real newlink = links( lat ) , oldlink , max ;
  *theta = theta_test_lin( lat , &max , ND ) ; 

  // initialise the draughtboard
  init_cb( LVOLUME ) ;

  printf( "[GF] Over-Relaxation parameter %f \n" , OrParam ) ;

  // iterations
  int iters = 0 ;
  while( iters < MAX_ITERS && fabs( *theta ) > ACC ) {

    oldlink = newlink ;

    // loop su2 indices
    int su2_index ;
    for( su2_index = 0 ; su2_index < NSU2SUBGROUPS ; su2_index++ ) {
      OR_iteration( lat , su2_index , OrParam , 0 , ND ) ;
    }

    newlink = links( lat ) ;

    printf( "%1.12f \n" , newlink ) ;

    // chroma condition is pretty shitty
    *theta = ( newlink - oldlink ) / newlink ;

    iters++ ;
  }

  // if not completed, grab file and redo

  // and print it out
  output_fixing_info( lat , *theta , iters ) ;

  free_cb( ) ;

  return iters ;
}

// little function for computing the spatial links for the Coulomb
// gauge fixing code
static double
slice_spatial_links( const struct site *__restrict lat ,
		     const int t )
{
  double sum = 0.0 ;
  int i ;
  #pragma omp parallel for private(i) reduction(+:sum)
  for( i = LCU*t ; i < LCU*(t+1) ; i++ ) {
    register double loc_sum = 0.0 ;
    int mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      loc_sum += (double)creal( trace( lat[i].O[mu] ) ) ;
    }
    sum = sum + (double)loc_sum ;
  }
  return sum / ( LCU * (ND-1) * NC ) ;
}

// Coulomb gauge fix -> need to think about restarting this
int
OrCoulomb( struct site *__restrict lat ,
	   double *theta ,
	   const int MAX_ITERS , 
	   const double ACC ,
	   const double OrParam )
{
  // initialise the draughtboarding
  init_cb( LCU ) ;

  printf( "[GF] Over-Relaxation parameter %f \n\n" , OrParam ) ;

  int t , iters = 0 ;
  for( t = 0 ; t < Latt.dims[ND-1] ; t++ ) {

    double newlink = slice_spatial_links( lat , t ) , oldlink ;
    int loc_iters = 0 ;
    *theta = 1.0 ;

    // iterations
    while( loc_iters < MAX_ITERS && fabs( *theta ) > ACC ) {

      oldlink = newlink ;

      // loop su2 indices
      int su2_index ;
      for( su2_index = 0 ; su2_index < NSU2SUBGROUPS ; su2_index++ ) {
	OR_iteration( lat , su2_index , OrParam , t , ND-1 ) ;
      }

      newlink = slice_spatial_links( lat , t ) ;

      // chroma condition is pretty shitty
      *theta = ( newlink - oldlink ) / newlink ;

      loc_iters++ ;
    }

    printf( "[GF] Slice :: %d {Stopped by convergence} \n[GF] Accuracy :: %1.5e"
	    " || Iterations :: %d\n[GF] Failures :: %d\n\n" , 
	    t , *theta , loc_iters , 0 ) ; 
    iters += loc_iters ;
  }

  // and print it out
  double splink , tlink ;
  all_links( lat , &splink , &tlink ) ;
  printf( "[GF] Tuning :: %f || Iterations :: %d ||\n[GF] Final Tlink :: %1.15f "
	  "|| Slink :: %1.15f \n[GF] Plaquette :: %1.15f \n" , 
	  Latt.gf_alpha , iters , tlink , splink , av_plaquette( lat ) ) ; 
  
  // memory frees
  free_cb( ) ;

  return iters ;
}
