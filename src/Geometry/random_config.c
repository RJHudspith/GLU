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
  GLU_complex *temp = malloc( NCNC * sizeof( GLU_complex ) ) ; 
  int i ; 
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    reunit2( gauge[i] ) ; 
  }

  free( temp ) ; 
  return ;
}

//lattice reunitarization for lattice links
void 
latt_reunitU( struct site *__restrict lat )
{
  int i ; 
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      reunit2( lat[i].O[mu] ) ; 
    }
  }
  return ;
}

// random gauge transform for our su(NC) matrix
void 
random_gtrans( struct site *__restrict lat )
{
  rng_init(  ) ; 
  
  printf( "\n[RNG] Performing a RANDOM gauge transformation \n" ) ;

  GLU_complex **gauge = malloc ( LVOLUME * sizeof( GLU_complex * ) ) ;
  int i ; 
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    gauge[i] = malloc ( NCNC * sizeof( GLU_complex ) ) ;
  }

  // openmp does not play nice with static arrays used in the WELL and MWC_1038
  for( i = 0 ; i < LVOLUME ; i++ ) { generate_NCxNC( gauge[i] ) ; }
 
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {    
    #if NC > 5 && ( ( defined HAVE_GSL ) || ( defined HAVE_LAPACKE_H ) )
    GLU_complex A[ NCNC ] ;
    Hermitian_proj( A , gauge[i] ) ;
    exponentiate( gauge[i] , A ) ;
    #else
    reunit2( gauge[i] ) ;
    #endif
  }

  gtransform( lat , ( const GLU_complex ** )gauge ) ; 

#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      reunit2( lat[i].O[mu] ) ; 
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
  int i ;
  // openmp does not play nice with RNG
  for( i = 0 ; i < LCU ; i ++  ) {
    generate_NCxNC( slice_gauge[i] ) ;
  } 	      
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i ++  ) {
    reunit2( slice_gauge[i] ) ;
  }
  return ;
}

// perform a random transform of the gauge links
void
random_transform( struct site *__restrict lat ,
		  GLU_complex *__restrict *__restrict gauge )
{
  int i ; 
  rng_init( ) ; 
  // openmp does not play nice with RNG
  for( i = 0 ; i < LVOLUME ; i++ ) {
    generate_NCxNC( gauge[i] ) ;
  }
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    reunit2( gauge[i] ) ;
  }
  gtransform( lat , ( const GLU_complex ** )gauge ) ; 
  return ;
}

// random gauge matrices
void 
random_suNC( struct site *__restrict lat )
{
  int i ;
  rng_init(  ) ; 
  // openmp does not play nice with RNG
  for( i = 0 ; i < LVOLUME ; i++ ) {
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      generate_NCxNC( lat[i].O[mu] ) ;
    }
  }
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {      
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      reunit2( lat[i].O[mu] ) ;
    }
  }
  return ;
}

/// reunit two slices
void
reunit_gauge_slices( GLU_complex **gauge1 , 
		     GLU_complex **gauge2 )
{
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    reunit2( gauge1[i] ) ;
    reunit2( gauge2[i] ) ;
  }
  return ;
}

// All ones 
void
trivial( struct site *__restrict lat )
{
  int i ; 
  printf( "\n[UNIT] Creating identity SU(%d) lattice fields \n" , NC ) ; 
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) { identity( lat[i].O[mu] ) ; } 
  }    
  return ;
}
