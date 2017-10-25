/*
    Copyright 2013-2016 Renwick James Hudspith

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

#include "par_rng.h"      // for generate_NCxNC()
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

//lattice reunitarization for the links
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
  initialise_par_rng( NULL ) ; 
 
  fprintf( stdout , "\n[RNG] Performing a RANDOM gauge transformation \n" ) ;

  GLU_complex **gauge = NULL ;
  
  if( GLU_malloc( (void**)&gauge , ALIGNMENT , LVOLUME * sizeof( GLU_complex* ) ) != 0 ) {
    fprintf( stdout , "[RANDOM] Allocation failed \n" ) ;
  }

  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_malloc( (void**)&gauge[i] , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) ;
    Sunitary_gen( gauge[i] , get_GLU_thread( ) ) ;
  }

  gtransform( lat , (const GLU_complex **)gauge ) ; 

#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    free( gauge[i] ) ;
  }
  free( gauge ) ;

  return ;
}

// randomly generate a new gauge transform slice
void
random_gtrans_slice( GLU_complex *__restrict *__restrict slice_gauge )
{
  initialise_par_rng( NULL ) ; 

  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i ++  ) {
    Sunitary_gen( slice_gauge[i] , get_GLU_thread( ) ) ;
  }
  
  return ;
}

// perform a random transform of the gauge links
void
random_transform( struct site *__restrict lat ,
		  GLU_complex *__restrict *__restrict gauge )
{
  initialise_par_rng( NULL ) ; 

  size_t i ; 
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    Sunitary_gen( gauge[i] , get_GLU_thread( ) ) ;
  }
  gtransform( lat , (const GLU_complex **)gauge ) ; 
  return ;
}

// random gauge matrices
void 
random_suNC( struct site *__restrict lat )
{
  initialise_par_rng( NULL ) ; 

  size_t i ;
  // openmp does not play nice with RNG
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    uint32_t thread = get_GLU_thread( ) ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      Sunitary_gen( lat[i].O[mu] , thread ) ;
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
