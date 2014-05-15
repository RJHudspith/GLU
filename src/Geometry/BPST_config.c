/*
    Copyright 2013 Renwick James Hudspith

    This file (BPST_config.c) is part of GLU.

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
   @file BPST_config.c
   @brief creates a BPST instanton configuration
 */
#include "Mainfile.h"
#include "geometry.h"

// numbers that seem to work rho = 4.0 (singular) rho = 1.0 (regular)
#define SINGULAR

#if ND == 4
static GLU_real
Asing( GLU_real y[ ND ] , const GLU_real rho , int mu )
{
  GLU_real sumr = 0.0 , sum = 0.0 ;
  int nu ;
  for( nu = 0 ; nu < ND ; nu++ ){
    if( nu != mu ){ sum += y[nu]*y[nu] ; }
    sumr += y[nu] * y[nu] + y[nu];
  }
  sum = sqrt(sum + rho*rho);
  sumr = sumr + rho*rho;
  
  return atan2( sum , sumr ) / sum; 
  //return atan( sum / sumr ) / sum; 
}
#endif

// create a BPST instanton lattice configuration
void
instanton_config( struct site *lat ) 
{
#if ND != 4
  printf( "Instanton solution not available for ND = %d \n" , ND ) ;
  return ;
#else
  // embed pauli matrices in SU(NC) matrices
  GLU_complex t[ 3 ][ NCNC ] ;  
  int mu ;
  for( mu = 0 ; mu < 3 ; mu++ ) { zero_mat( t[mu] ) ; }
  t[0][1] = 1. ; t[0][NC] = 1.0 ;
  t[1][1] = -I ; t[1][NC] = I ;
  t[2][0] = 1.0 ; t[2][NC+1] = -1.0 ;

  const GLU_real rho = 4.0 ; // instanton size
  const GLU_real sign = 1 ; // instanton (-1) or anti-instanton (1)
  GLU_real c[ ND ] ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // put it just off the centre
    c[ mu ] = ( Latt.dims[mu]/2 - 0.5 ) ;
  }
  //const GLU_real c[ ND ] = { 3.5 , 3.5 , 3.5 , 3.5 } ; // centre
  //const GLU_real c[ ND ] = { 3 , 3 , 3 , 3 } ; // centre

  int i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {

    // compute lattice position vector ...
    int x[ ND ] ;
    get_mom_2piBZ( x , i , ND ) ;

    // compute the coordinate "y"
    GLU_real y[ ND ] ;
    int nu ;
    for( nu = 0 ; nu < ND ; nu++ ) {
      y[ nu ] = x[ nu ] - c[ nu ] ;
    }
    
    // compute "b"
    GLU_real b[ ND ][ 3 ] ;
    b[0][0] = -sign *y[3]; b[0][1] =  y[2];       b[0][2] = -y[1];
    b[1][0] = -y[2];       b[1][1] = -sign *y[3]; b[1][2] =  y[0];
    b[2][0] =  y[1];       b[2][1] = -y[0];       b[2][2] = -sign *y[3];
    b[3][0] =  sign *y[0]; b[3][1] =  sign *y[1]; b[3][2] =  sign *y[2]; 
 
    // perform the exponentiation, loop nu
    for( nu = 0 ; nu < ND ; nu++ ) {
      GLU_complex temp[ NCNC ] = {} ;
      int j ;
      for( j = 0 ; j < 3 ; j++ ) {
	int k ;
	for( k = 0 ; k < NCNC ; k++ ) {
	  temp[ k ] += b[ nu ][ j ] * t[ j ][ k ] ; 
	}
      }
      // multiply by asing
      #ifdef SINGULAR
      const GLU_real c = Asing( y , 0 , nu ) - Asing( y , rho , nu ) ;
      #else
      const GLU_real c = Asing( y , rho , nu ) ;
      #endif
      for( j = 0 ; j < NCNC ; j++ ) {
	temp[ j ] *= c ; 
      }
      exponentiate( lat[i].O[nu] , temp ) ;
    }
  }
  return ;
#endif
}
