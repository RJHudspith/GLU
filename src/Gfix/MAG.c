/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (MAG.c) is part of GLU.

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
   @file MAG.c
   @brief Maximal (temporal) Axial Gauge fixing

  Here I have included the Maximal Axial Gauge fixing and
  the \f$ A^{\mu} = 0 \f$ partial axial gauge fixing.

  MAG performs a lattice-wide reunitarisation step first.
 */
#include "Mainfile.h"

#include "geometry.h"      // get site 
#include "givens.h"        // residual fixing uses this
#include "gtrans.h"        // lattice-wide gauge transformations
#include "gramschmidt.h"   // orthogonalisation
#include "random_config.h" // lattice reunitarisation and stuff

//calculate gauge transform matrices along one line in one direction
static void 
gauge_set( GLU_complex **gauge ,
	   const struct site *lat , 
	   const size_t mu , 
	   const size_t i )
{  
  size_t left = i ; 
  size_t next = lat[i].neighbor[mu] ; 
  // set origin to identity
  // I guess I could just equate this with lat to save on a multiply
  identity( gauge[i] ) ; 
  //loop along the lattice
  size_t j ; 
  for( j = 0 ; j < ( Latt.dims[ mu ]-1 ) ; j++ ) {
    multab_suNC( gauge[ next ] , gauge[left] , lat[left].O[mu] ) ; 
    left = next ; 
    next = lat[ left ].neighbor[mu] ; 
  }
  return ;
}

//perform the maximal axial gauge fixing
static int 
MAG_fix( struct site *lat )
{
  // allocate gauge field
  GLU_complex **gauge = NULL ;
  size_t i ;
  int flag = GLU_SUCCESS ;
  
  // allocate new gauge field
  if( GLU_malloc( (void**)&gauge , 16 , LVOLUME * sizeof( GLU_complex* ) ) != 0 ) {
    fprintf( stderr , "[GF] MAG failed to allocate temporary gauge\n" ) ;
    flag = GLU_FAILURE ;
    goto end ;
  }
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    gauge[i] = ( GLU_complex* )malloc( NCNC * sizeof( GLU_complex ) ) ; 
    identity( gauge[i] ) ;
  }

  // set the gauge transformation matrices
  size_t subvol = 1 , mu ;
  for( mu = 0 ; mu < ND ; mu ++ ) {
    #pragma omp parallel for private(i)
    for( i =  0 ; i < subvol ; i++ ) {
      gauge_set( gauge , lat , mu , i ) ;
    } 
    gtransform( lat , ( const GLU_complex ** )gauge ) ;
    subvol *= Latt.dims[ mu ] ;
  }

 end :

  // free the gauge transformation matrices
  if( gauge != NULL ) {
#pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      free( gauge[i] ) ;   
    }
    free( gauge ) ;
  }
  
  return flag ;
}

// simplistic axial gauge calculation, A_{DIR} = 0 as best we can
int
axial_gauge( struct site *lat ,
	     const size_t DIR )
{
  latt_reunitU( lat ) ; 
  
  GLU_complex **gauge ;
  size_t i , subvolume = 1 ;
  int flag = GLU_SUCCESS ;
  
  // allocate new gauge field
  if( GLU_malloc( (void**)&gauge , 16 , LVOLUME * sizeof( GLU_complex* ) ) != 0 ) {
    fprintf( stderr , "[GF] Axial failed to allocate temporary gauge\n" ) ;
    flag = GLU_FAILURE ;
    goto end ;
  }
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    gauge[i] = ( GLU_complex* )malloc( NCNC * sizeof( GLU_complex ) ) ; 
    identity( gauge[i] ) ;
  }
  
  for( i =  0 ; i < ND ; i++ ) {
    subvolume *= (i!=DIR) ? Latt.dims[i] : 1 ;
  }
#pragma omp parallel for private(i)
  for( i =  0 ; i < subvolume ; i++ ) {// loop the ND-1, subvolume
    // x is the ND dimensional vector describing the position
    int x[ ND ] ;
    get_mom_2piBZ( x , i , DIR ) ;
    const int k = gen_site( x ) ;
    gauge_set( gauge , lat , DIR , k ) ; 
  }
  gtransform( lat , ( const GLU_complex ** )gauge ) ;

 end :
  
  // free the gauge transformation matrices
  if( gauge != NULL ) {
#pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      free( gauge[i] ) ;   
    }
    free( gauge ) ;
  }

  latt_reunitU( lat ) ; 
  
  return flag ;
}

// wrapper for the mag-fixing bit ...
int
mag( struct site *lat )
{
  //loop lattice reuinitarizing everything
  latt_reunitU( lat ) ; 
  int flag = MAG_fix( lat ) ; 
  latt_reunitU( lat ) ; 
  return flag ;
}

// residual fix
int
residual_fix( struct site *lat )
{
  // temporal gauge transformation matrices
  GLU_complex **gauge = NULL ;
  size_t i ;
  int flag = GLU_SUCCESS ;
  
  // set these
  if( GLU_malloc( (void**)&gauge , 16 , Latt.dims[ND-1] * sizeof( GLU_complex* ) ) != 0 ) {
    fprintf( stderr , "[RESFIX] temporary gauge allocation failure\n" ) ;
    flag = GLU_FAILURE ;
    goto end ;
  }

  for( i = 0 ; i < Latt.dims[ND-1] ; i++ ) {
    if( GLU_malloc( (void**)&gauge[i] , 16 , NCNC * sizeof( GLU_complex ) ) != 0 ) {
      fprintf( stderr , "[RESFIX] temporary gauge allocation failure\n" ) ;
      flag = GLU_FAILURE ;
      goto end ;
    }
    if( i == 0 ) identity( gauge[0] ) ;
  } 

  const GLU_real one_LCU = 1.0 / LCU ; 
  size_t t ;
  for( t = 0 ; t < Latt.dims[ND-1]-1 ; t++ ) {
    GLU_complex sum[ NCNC ] GLUalign ;
    zero_mat( sum ) ;
    for( i = 0 ; i < LCU ; i++ ) {
      a_plus_b( sum , lat[ i + LCU*t ].O[ND-1] ) ;
    }
    // turn it into the average ... not needed?
    for( i = 0 ; i < NCNC ; i++ ) { sum[i] *= one_LCU ; }

    // project the sum of the temporal links back to SU(NC) 
    // two ways of doing this I can think of
    #ifdef HERM_PROJ
    // 1. Hermitian projection and exact exponentiation
    GLU_complex A[ NCNC ] GLUalign ;
    Hermitian_proj( A , sum ) ;
    exponentiate( sum , A ) ;
    #else
    // 2. Simple reunitarisation
    givens_reunit( sum ) ;
    #endif

    // and multiply the next gauge transformations with the 
    // previous and this one!
    multab_suNC( gauge[t+1] , gauge[t] , sum ) ;
  } 
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex temp[ NCNC ] GLUalign ;
    const size_t t = (int)( i / LCU ) ;
    const size_t tup = ( t == Latt.dims[ND-1]-1 ) ? 0 : t+1 ;
    multab_dag_suNC( temp , lat[i].O[ND-1] , gauge[tup] ) ;
    multab_suNC( lat[i].O[ND-1] , gauge[t] , temp ) ; 
    #ifdef SINGLE_PREC
    reunit2( lat[i].O[ND-1] ) ;
    #endif
    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu ++ ) {
      multab_dag_suNC( temp , lat[i].O[mu] , gauge[t] ) ;
      multab_suNC( lat[i].O[mu] , gauge[t] , temp ) ; 
      #ifdef SINGLE_PREC
      reunit2( lat[i].O[mu] ) ;
      #endif
    }
  }

 end :

  // free the temporary gauges
  if( gauge != NULL ) {
    for( i = 0 ; i < Latt.dims[ND-1] ; i++ ) {
      free( gauge[i] ) ;
    }
    free( gauge ) ;
  }
  
  return flag ;
}
