/*
    Copyright 2013 Renwick James Hudspith

    This file (Qsusc.c) is part of GLU.

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
   @file Qsusc.c
   @brief computation of the topological susceptibility correlator
   @warning calls the smearing wrapper, which overwrites the links
 */

#include "Mainfile.h"

#include "clover.h"       // computation of the topological charge
#include "cut_output.h"   // automatic formatting of our output file
#include "cut_routines.h" // momentum cuts for config-space vector
#include "GLU_bswap.h"    // byte swapping
#include "SM_wrap.h"      // in case we wish to do smearing

// computation of the correlator
// C(r) = < q( x ) q( y ) >  -> r = ( x - y ), r^2 < CUTINFO.max_mom
// over all indices, q(x) is the topological charge at site x
void
compute_Qsusc( struct site *__restrict lat ,
	       const struct cut_info CUTINFO ,
	       const struct sm_info SMINFO )
{
  // normalisations
  const double NORM = -0.001583143494411527678811 ; // 1.0/(64*Pi*Pi)
  const double NORMSQ = NORM * NORM ;

  // do some smearing ...
  SM_wrap_struct( lat , SMINFO ) ;

  // info
  check_psq( CUTINFO ) ;

  // compute the ratios of the dimensions in terms of the smallest
  simorb_ratios( ND ) ;

  // if we already have a file, read it
  int size[1] = {} ;
  struct veclist *list = compute_veclist( size , CUTINFO , ND , GLU_TRUE ) ;

  // set up the outputs
  const char *str = output_str_struct( CUTINFO ) ;  
  FILE *Ap = fopen( str , "wb" ) ; 

  // write the list, in cut_outputs.c
  write_mom_veclist( Ap , size , list , ND ) ;

  // set up the matrix-valued array of the topological charge
  GLU_complex **qtop = malloc( LVOLUME * sizeof( GLU_complex* ) ) ;

  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    qtop[i] = (GLU_complex*)malloc( NCNC * sizeof( GLU_complex ) ) ;
  }

  // precompute all of charge densities q(x)
  compute_Gmunu_array( qtop , lat ) ;

  // allocate the results
  double *qcorr = malloc( size[0] * sizeof( double ) ) ; 
  double *trtr  = malloc( size[0] * sizeof( double ) ) ; 

  // loop the possible rsqs
  #pragma omp parallel for private(i)
  for( i = 0 ; i < size[0] ; i++ ) {
    // some storage for the traces
    GLU_complex tr , tr1 , tr2 ;
    register double sum = 0.0 , sumtrtr = 0.0 ;
    // loop the volume, performing a translationally invariant sum
    int j , sep ;
    for( j = 0 ; j < LVOLUME ; j++ ) {
      const int tmp = list[i].idx + j ;
      sep = tmp < LVOLUME ? tmp : tmp - LVOLUME ;
      //trace of the products
      trace_ab_dag( &tr , qtop[ sep ] , qtop[ j ] ) ;
      sum += creal( tr ) ;
      // compute the trace-trace portion
      tr1 = trace( qtop[ sep ] ) ;
      tr2 = trace( qtop[ j ] ) ;
      sumtrtr += creal( tr1 ) * creal( tr2 ) + 
	         cimag( tr1 ) * cimag( tr2 ) ;
    }
    qcorr[i] = sum * NORMSQ ;
    trtr[i]  = sumtrtr * NORMSQ ;
  }

  // write out the result of all that work
  write_g2g3_to_list( Ap , qcorr , trtr , size ) ;

  // close the file and its name
  fclose( Ap ) ;
  free( (void*)str ) ;

  // free up all of the allocations
  free( qcorr ) ;  
  free( trtr ) ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    free( qtop[i] ) ;
  }
  free( qtop ) ;

  // free the list
  free( list ) ;

  return ;
}
