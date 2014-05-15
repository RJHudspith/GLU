/*
    Copyright 2013 Renwick James Hudspith

    This file (config_gluons.c) is part of GLU.

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
   @file config_gluons.c
   @brief code to compute the configuration space gluonic correlation functions
 */

#include "Mainfile.h"

#include "cut_output.h"   // output file
#include "cut_routines.h" // momentum cutting?
#include "SM_wrap.h"      // do we want some smearing?

/**
   @param LT
   @brief length of the temporal direction
 */
#define LT Latt.dims[ND-1]

// compute the log fields, overwrite the link matrices!
static int
Amu_fields( A , def )
     struct site *__restrict A ;
     const lie_field_def def ;
{
  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    struct spt_site temp ;
    int mu ;
    switch( def )
      {
      case LINEAR_DEF :
	for( mu = 0 ; mu < ND ; mu++ ) 
	  Hermitian_proj( temp.O[mu] , A[i].O[mu] ) ;
	break ;
      case LOG_DEF :
	for( mu = 0 ; mu < ND ; mu++ ) 
	  exact_log_slow( temp.O[mu] , A[i].O[mu] ) ;
	break ;
      }
    memcpy( &A[i] , &temp , sizeof( struct spt_site ) ) ;
  }
  return GLU_SUCCESS ;
}

/////////// spatial-spatial point sources over the whole lattice  ////////
static void
spatial_pointsource( A , gs , gsnorm )
     const struct site *__restrict A ;
     double *__restrict gs ;
     const double gsnorm ;
{
  // loop time 
  int t ;
#pragma omp parallel for private(t)
  for( t = 0 ; t < LT/2+1 ; t++ ) { // is exactly symmetric around LT/2
    // local trace
    register double loc_tr = 0.0 ;
    // temporary trace
    GLU_real tr ;
    // loop separations tau
    int tau , i , j ;
    for( tau = 0 ; tau < LT ; tau++ ) {
      const int idx1 = LCU * tau ;
      const int idx2 = LCU * ( ( tau + t ) % LT ) ;
      // loop spatial hypercube index i
      for( i = 0 ; i < LCU ; i++ ) {
	const int x = i + idx1 ;
	for( j = 0 ; j < LCU ; j++ ) {
	  const int y = j + idx2 ;
	  #if ND == 4
	  trace_ab_herm( &tr , A[x].O[0] , A[y].O[0] ) ;
	  loc_tr += (double)tr ;
	  trace_ab_herm( &tr , A[x].O[1] , A[y].O[1] ) ;
	  loc_tr += (double)tr ;
	  trace_ab_herm( &tr , A[x].O[2] , A[y].O[2] ) ;
	  loc_tr += (double)tr ;
	  //trace_ab_herm( &tr , A[x].O[3] , A[y].O[3] ) ;
	  //loc_tr += (double)tr ;
	  #else
	  int mu ;
	  for( mu = 0 ; mu < ND-1 ; mu++ ) {
	    trace_ab_herm( &tr , A[x].O[mu] , A[y].O[mu] ) ;
	    loc_tr += (double)tr ;
	  }
	  #endif
	}
      }
    }
    gs[t] = loc_tr * gsnorm / LT ;
  }
#pragma omp parallel for private(t)
  for( t = LT/2+1 ; t < LT ; t++ ) {
    // all to all is symmetric around LT/2
    gs[ t ] = gs[ (LT)-t ] ;
  }
  return ;
}

// computes the configuration space propagators
int 
cuts_struct_configspace( struct site *__restrict A ,
			 const struct cut_info CUTINFO ,
			 const struct sm_info SMINFO )
{
  // smeared gluon field to extract the ground state?
  SM_wrap_struct( A , SMINFO ) ;

  // take the log
  if( Amu_fields( A , CUTINFO.definition ) != GLU_SUCCESS ) {
    printf( "[CUTS] Something wrong in Logging of the fields, check it out \n" ) ;
    return GLU_FAILURE ;
  } 

  printf( "\n[CUTS] Computing the configuration-space gluon propagators.\n" ) ;

  // set up the outputs
  char *str = output_str_struct( CUTINFO ) ;  
  FILE *Ap = fopen( str , "wb" ) ; 

  // compute the temporal and spatial correlators
  double *gsp = malloc( Latt.dims[ ND-1 ] * sizeof( double ) ) ;

  // normalisations
  const double spat_norm = 1.0 / ( ( ND-1 ) * LCU ) ;

  // choices, choices. I have now gone for wall point and the wall-wall
  spatial_pointsource( A , gsp , spat_norm ) ;

  int lt[ 1 ] = { LT } ;

  // write out the timeslice list ...
  write_tslice_list( Ap , lt ) ;

  // and write the props ....
  write_g2_to_list( Ap , gsp , lt ) ; 

  free( gsp ) ;

  return GLU_SUCCESS ;
}

#undef LT // clean it up?
