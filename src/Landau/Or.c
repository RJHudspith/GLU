/*
    Copyright 2013-2016 Renwick James Hudspith

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

   TODO :: is this correct? Need to think - J
 */

#include "Mainfile.h"       // general includes

#include "draughtboard.h"   // draughtboarding
#include "gftests.h"        // theta test
#include "gtrans.h"         // gauge transformations
#include "plaqs_links.h"    // plaquettes
#include "random_config.h"  // latt reunitisation
#include "SU2_rotate.h"     // su(2) rotations

#if 0
// NC generic givens rotations
static void
OrRotation( GLU_complex *__restrict s0 , 
	    GLU_complex *__restrict s1 ,
	    const GLU_complex U[ NCNC ] , 
	    const double OrParam ,
	    const size_t su2_index )
{
  *s0 = U[ Latt.su2_data[ su2_index ].idx_a ] 
    + conj( U[ Latt.su2_data[ su2_index ].idx_d ] ) ;
  *s1 = U[ Latt.su2_data[ su2_index ].idx_b ] 
    - conj( U[ Latt.su2_data[ su2_index ].idx_c ] ) ;

  // I use the MILC OR step here
  register const GLU_real asq = cimag(*s0)*cimag(*s0) 
    + creal(*s1)*creal(*s1) + cimag(*s1)*cimag(*s1) ;
  register const GLU_real a0sq = creal(*s0)*creal(*s0) ;
  register const GLU_real x = ( OrParam * a0sq + asq ) / ( a0sq + asq ) ;
  register const GLU_real r = sqrt( a0sq + x*x*asq ) ;
  register const GLU_real xdr = x/r ;
  *s0 = creal( *s0 ) / r + I * cimag( *s0 ) * xdr ;
  *s1 *= xdr ;

  return ;
}

// a single iteration of the overrelaxed gauge fixing
static void
OR_single( struct site *__restrict lat ,
	   const size_t su2_index ,
	   const double OrParam ,
	   const size_t i ,
	   const size_t DIMS )
{
  GLU_complex L[ NCNC ] GLUalign ;
  zero_mat( L ) ;
  size_t mu , j , k ;
  // loop directions summing into L
  for( mu = 0 ; mu < DIMS ; mu++ ) {
    // compute U(x+\mu/2) + U^{dagger}(x-\mu/2)
    const size_t back = lat[i].back[mu] ;
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
  OrRotation( &s0 , &s1 , L , OrParam , su2_index ) ;

  // gauge rotate
  for( mu = 0 ; mu < ND ; mu++ ) {
    const size_t back = lat[i].back[mu] ;
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
	      const struct draughtboard db ,
	      const size_t su2_index ,
	      const double OrParam ,
	      const size_t t ,
	      const size_t DIMS )
{
  // perform an overrelaxation step
  size_t i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < db.Nred ; i++ ) {  
    OR_single( lat , su2_index  , OrParam , 
	       db.red[i] + LCU * t , DIMS ) ;
  }
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < db.Nblack ; i++ ) {    
    OR_single( lat , su2_index , OrParam , 
	       db.black[i] + LCU * t , DIMS ) ;
  }
  return ;
}

// output the data, pass lat for the plaquette
static void
output_fixing_info( struct site *__restrict lat ,
		    const double theta ,
		    const size_t iters )
{
  // reunitarise just to limit the damage from round-off
  latt_reunitU( lat ) ;

  ////////// Print out the Gauge Fixing information /////////////
  fprintf( stdout , "[GF] Plaquette :: %1.15f \n[GF] Accuracy :: %1.4e\n" , 
	   av_plaquette( lat ) , theta ) ;
  GLU_real tr ;
  const double link = indivlinks( lat , &tr ) ;
  fprintf( stdout , "[GF] Iters :: %zu\n[GF] Link trace :: %1.15f ||"
	   "Maximum :: %1.15f\n" , iters , link , tr / NC ) ; 
  double lin , log ;
  const_time( lat , &lin , &log ) ; 
  fprintf( stdout , "[GF] Temporal constance || Lin %e || Log %e \n" , 
	   lin , log ) ;
  gauge_functional( lat ) ;
  fprintf( stdout , "[GF] Functional :: %1.15f\n" , gauge_functional( lat ) ) ;
  ///////////////////////////////////////////////////////////////
  return ;
}

// gauge fix
size_t
OrLandau( struct site *__restrict lat ,
	  double *theta ,
	  const size_t MAX_ITERS , 
	  const double ACC ,
	  const double OrParam )
{
  GLU_real newlink = links( lat ) , oldlink , max ;
  *theta = theta_test_lin( lat , &max , ND ) ; 

  // initialise the draughtboard
  struct draughtboard db ;
  init_cb( &db , LVOLUME , ND ) ;

  fprintf( stdout , "[GF] Over-Relaxation parameter %f \n" , OrParam ) ;

  // iterations
  size_t iters = 0 ;
  while( iters < MAX_ITERS && fabs( *theta ) > ACC ) {

    oldlink = newlink ;

    // loop su2 indices
    size_t su2_index ;
    for( su2_index = 0 ; su2_index < NSU2SUBGROUPS ; su2_index++ ) {
      OR_iteration( lat , db , su2_index , OrParam , 0 , ND ) ;
    }

    newlink = links( lat ) ;

    #ifdef verbose
    fprintf( stdout , "%1.12f \n" , newlink ) ;
    #endif

    // chroma condition is pretty shitty
    *theta = ( newlink - oldlink ) / newlink ;

    iters++ ;
  }

  // and print it out
  output_fixing_info( lat , *theta , iters ) ;

  free_cb( &db ) ;

  return iters ;
}

// little function for computing the spatial links for the Coulomb
// gauge fixing code
static double
slice_spatial_links( const struct site *__restrict lat ,
		     const size_t t )
{
  double sum = 0.0 ;
  size_t i ;
  #pragma omp parallel for private(i) reduction(+:sum)
  for( i = LCU*t ; i < LCU*(t+1) ; i++ ) {
    register double loc_sum = 0.0 ;
    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      loc_sum += (double)creal( trace( lat[i].O[mu] ) ) ;
    }
    sum = sum + (double)loc_sum ;
  }
  return sum / ( LCU * (ND-1) * NC ) ;
}

// Coulomb gauge fix -> need to think about restarting this
size_t
OrCoulomb( struct site *__restrict lat ,
	   double *theta ,
	   const size_t MAX_ITERS , 
	   const double ACC ,
	   const double OrParam )
{
  // initialise the draughtboarding
  struct draughtboard db ;
  init_cb( &db , LCU , ND-1 ) ;

  fprintf( stdout , "[GF] Over-Relaxation parameter %f \n\n" , OrParam ) ;

  size_t t , iters = 0 ;
  for( t = 0 ; t < Latt.dims[ND-1] ; t++ ) {

    double newlink = slice_spatial_links( lat , t ) , oldlink ;
    size_t loc_iters = 0 ;
    *theta = 1.0 ;

    // iterations
    while( loc_iters < MAX_ITERS && fabs( *theta ) > ACC ) {

      oldlink = newlink ;

      // loop su2 indices
      size_t su2_index ;
      for( su2_index = 0 ; su2_index < NSU2SUBGROUPS ; su2_index++ ) {
	OR_iteration( lat , db , su2_index , OrParam , t , ND-1 ) ;
      }

      newlink = slice_spatial_links( lat , t ) ;

      // chroma condition is pretty shitty
      *theta = ( newlink - oldlink ) / newlink ;

      loc_iters++ ;
    }

    fprintf( stdout , "[GF] Slice :: %zu {Stopped by convergence} \n"
	     "[GF] Accuracy :: %1.5e || Iterations :: %zu\n"
	     "[GF] Failures :: %d\n\n" , t , *theta , loc_iters , 0 ) ; 
    iters += loc_iters ;
  }

  // and print it out
  double splink , tlink ;
  all_links( lat , &splink , &tlink ) ;
  fprintf( stdout , "[GF] Tuning :: %f || Iterations :: %zu ||\n"
	   "[GF] Final Tlink :: %1.15f || Slink :: %1.15f \n"
	   "[GF] Plaquette :: %1.15f \n" , Latt.gf_alpha , iters , 
	   tlink , splink , av_plaquette( lat ) ) ; 
  // memory frees
  free_cb( &db ) ;

  return iters ;
}

#endif
