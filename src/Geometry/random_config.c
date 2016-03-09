/*
    Copyright 2013 Renwick James Hudspith

    This file (random_config.c) is part of GLU.

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
   @file random_config.c
   @brief create a few either random or trivial configurations or reunitarise some of our matrices
 */

#include "Mainfile.h"

#include "GLU_rng.h"      // for generate_NCxNC()
#include "gramschmidt.h"  // orthogonalisation
#include "gtrans.h"       // gauge transformations

//lattice reunitarization for gauge field
void 
latt_reunitg( GLU_complex *__restrict *__restrict gauge )
{
  size_t i ; 
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    gram_reunit( gauge[i] ) ; 
  }
  return ;
}

//lattice reunitarization for lattice links
void 
latt_reunitU( struct site *__restrict lat )
{
  size_t i ; 
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      gram_reunit( lat[i].O[mu] ) ; 
    }
  }
  return ;
}

// random gauge transform for our su(NC) matrix
void 
random_gtrans( struct site *__restrict lat )
{
  rng_init(  ) ; 
  
  fprintf( stdout , "\n[RNG] Performing a RANDOM gauge transformation \n" ) ;

  GLU_complex **gauge = NULL ;
  
  GLU_malloc( (void**)&gauge , ALIGNMENT , LVOLUME * sizeof( GLU_complex* ) ) ;

  size_t i ; 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_malloc( (void**)&gauge[i] , ALIGNMENT , LVOLUME * sizeof( GLU_complex ) ) ;
    // openmp does not play nice with static arrays in the WELL and MWC_1038
    generate_NCxNC( gauge[i] ) ;
    #if NC > 5 && ( ( defined HAVE_GSL ) || ( defined HAVE_LAPACKE_H ) )
    GLU_complex A[ NCNC ] ;
    Hermitian_proj( A , gauge[i] ) ;
    exponentiate( gauge[i] , A ) ;
    #else
    gram_reunit( gauge[i] ) ;
    #endif
  }
  gtransform( lat , (const GLU_complex **)gauge ) ; 
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      gram_reunit( lat[i].O[mu] ) ; 
    }
   free( gauge[i] ) ;
  }
  free( gauge ) ;

  return ;
}

// randomly gauge transform a slice
void
random_gtrans_slice( GLU_complex *__restrict *__restrict slice_gauge )
{
  rng_init( ) ; 
  size_t i ;
  // openmp does not play nice with RNG
  for( i = 0 ; i < LCU ; i ++  ) {
    generate_NCxNC( slice_gauge[i] ) ;
  } 	      
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i ++  ) {
    gram_reunit( slice_gauge[i] ) ;
  }
  return ;
}

// perform a random transform of the gauge links
void
random_transform( struct site *__restrict lat ,
		  GLU_complex *__restrict *__restrict gauge )
{
  size_t i ; 
  rng_init( ) ; 
  // openmp does not play nice with RNG
  for( i = 0 ; i < LVOLUME ; i++ ) {
    generate_NCxNC( gauge[i] ) ;
  }
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    gram_reunit( gauge[i] ) ;
  }
  gtransform( lat , (const GLU_complex **)gauge ) ; 
  return ;
}

// random gauge matrices
void 
random_suNC( struct site *__restrict lat )
{
  size_t i ;
  rng_init(  ) ; 
  // openmp does not play nice with RNG
  for( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      generate_NCxNC( lat[i].O[mu] ) ;
    }
  }
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {      
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      gram_reunit( lat[i].O[mu] ) ;
    }
  }
  return ;
}

/// reunit two slices
void
reunit_gauge_slices( GLU_complex **gauge1 , 
		     GLU_complex **gauge2 )
{
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    gram_reunit( gauge1[i] ) ;
    gram_reunit( gauge2[i] ) ;
  }
  return ;
}

// All ones 
void
trivial( struct site *__restrict lat )
{
  size_t i ; 
  fprintf( stdout , "\n[UNIT] Creating identity SU(%d) lattice fields \n" , 
	   NC ) ; 
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) { identity( lat[i].O[mu] ) ; } 
  }    
  return ;
}
