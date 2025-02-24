/*
Copyright 2013-2025 Renwick James Hudspith

    This file (CG.c) is part of GLU.

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
   @file CG.c
   @brief common code for the the CG routines
 */
#include "Mainfile.h"    // general definitions

#include "GLU_splines.h" // GLUbic spline interpolation code
#include "GLU_sums.h"    // round-off resistant summations
#include "gramschmidt.h" // for reunitarisation
#include "gtrans.h"      // gauge transformations

// this is the point where the gram-schmidt loses out
// to the n-ape det-rescaled projection
#if NC > 18
#include "taylor_logs.h"
#endif

// for the polyak-ribiere
double
PRfmax( const double a , const double b )
{
  return ( b < a ? a : b ) ;
}

// this is the same between Landau and Coulomb so I put it here
// I know that this probably cache misses really badly, the difficulty is
// that fftw wants the "in" matrix to be ordered this way
void
set_gauge_matrix( GLU_complex *gx ,
		  const GLU_complex **in ,
		  const double alpha ,
		  const size_t i )
{
#ifdef exp_exact
  GLU_complex short_gauge[ HERMSIZE ] ;
  #if NC==3
  short_gauge[ 0 ] = alpha * cimag( in[ 0 ][ i ] ) ;
  short_gauge[ 1 ] = -I * alpha * in[ 1 ][ i ] ;
  short_gauge[ 2 ] = -I * alpha * in[ 2 ][ i ] ;
  short_gauge[ 3 ] = -alpha * creal( in[ 0 ][ i ] ) ;
  short_gauge[ 4 ] = -I * alpha * in[ 3 ][ i ] ;
  #elif NC==2
  short_gauge[ 0 ] = alpha * cimag( in[ 0 ][ i ] ) ;
  short_gauge[ 1 ] = -I * alpha * in[ 1 ][ i ] ;
  #else 
  // this makes it hermitian, which is what the exponential expects
  int mu ;
  for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
    short_gauge[mu] = -I * alpha * in[ mu ][ i ] ;
  }
  #endif
  exponentiate_short( gx , short_gauge ) ;
#else
  #if NC==3
  gx[0] = 1.0 + I * alpha * cimag( in[0][i] ) ; 
  gx[1] = alpha * in[1][i] ;  
  gx[3] = -conj( gx[1] ) ; 
  gx[4] = 1.0 - I * alpha * creal( in[0][i] ) ; 
  gx[6] = -alpha * conj( in[2][i] ) ;  
  gx[7] = -alpha * conj( in[3][i] ) ;
  #elif NC==2
  gx[0] = 1.0 + I * alpha * cimag( in[0][i] ) ; 
  gx[1] = alpha * in[1][i] ; 
  gx[2] = -conj( gx[1] ) ; 
  #else
  GLU_complex short_gauge[ HERMSIZE ] ;
  // this makes it antihermitian
  size_t mu ;
  for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
    short_gauge[mu] = alpha * in[ mu ][ i ] ;
  }
  rebuild_antihermitian( gx , short_gauge ) ;
  add_constant( gx , 1.0 ) ; 
  #endif
  #if NC > 18
  nape_reunit( gx ) ; 
  #else
  gram_reunit( gx ) ; 
  #endif
#endif
  return ;
}
