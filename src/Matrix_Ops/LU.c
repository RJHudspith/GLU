/*
Copyright 2013-2025 Renwick James Hudspith

    This file (LU.c) is part of GLU.

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
   @file LU.c
   @brief LU decompositions
 */
#include "Mainfile.h"

// we have better determinant routines for small matrices
#if !( defined HAVE_IMMINTRIN_H ) || ( defined SINGLE_PREC ) && (NC>3)

#include "LU.h"        // alphabetising

//  the column-pivoted LU decomposition determinant
//  does not save L, just need the diagonal of U as determinant is product
//  of these elements
double complex
LU_det( const size_t N , const GLU_complex U[ N*N ] )
{
  double complex a[ N*N ] ;
  memcpy( a , U , N*N * sizeof( GLU_complex ) ) ;
  return LU_det_overwrite( N , a ) ;
}

// same as above, overwrites U
double complex
LU_det_overwrite( const size_t N , GLU_complex U[ N*N ] )
{
  size_t i , j , l , piv , perms = 0 ;
  double complex dt , determinant = 1. ;

  double attempt , best ; 
  for( i = 0 ; i < N-1 ; i++ ) {
    // swap rows s.t a[i] has largest pivot number first
    best = creal( U[i*(N+1)] ) * creal( U[i*(N+1)] ) 
         + cimag( U[i*(N+1)] ) * cimag( U[i*(N+1)] ) ;
    piv = i ;
    // again only care about the pivots below i
    for( j = i+1 ; j < N ; j++ ) {
      attempt = creal( U[i+j*N] ) * creal( U[i+j*N] ) 
	      + cimag( U[i+j*N] ) * cimag( U[i+j*N] ) ;
      if( attempt > best ) { 
	piv = j ; 
	best = attempt ; 
      }
    }
    if( best == 0.0 ) { 
      fprintf( stderr , "[DETERMINANT] LU  Singular Matrix!!!\n" ) ;
      return 0.0 ;
    }
    if( piv != i ) {
      // unlike what I am told, I physically swap rows
      // this is measured to be faster than saving the permutations
      // which I find quite weird, must be a caching thing
      for( l = i ; l < N ; l++ ) {
	dt         = U[l+i*N] ;
	U[l+i*N]   = U[l+piv*N] ;
	U[l+piv*N] = dt ;
      }
      perms++ ;
    }
    // perform gaussian elimination
    dt = 1.0 / U[ i*(N+1) ] ;
    double complex *pA = U + i*N ;
    for( j = N-1 ; j > i ; j-- ) { // go up in other column
      register double complex fac1 = U[ i + j*N ] * dt ; 
      // go along the row performing the subtraction, there is no point in
      // subtracting elements where we have determined the best pivot, just the
      // columns to the right of the pivot
      for( l = i + 1 ; l < N ; l++ ) {
	U[ l + j*N ] -= creal( fac1 ) * creal( pA[l] ) - cimag( fac1 ) * cimag( pA[l] ) 
	        + I * ( creal( fac1 ) * cimag( pA[l] ) + cimag( fac1 ) * creal( pA[l] ) ) ;
      }
    }
    determinant *= U[ i*(N+1) ] ;
  }
  determinant *= U[ N*N-1 ] ;
  return perms&1 ? -determinant : determinant ;
}

#endif
