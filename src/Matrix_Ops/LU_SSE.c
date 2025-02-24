/*
Copyright 2013-2025 Renwick James Hudspith

    This file (LU_SSE.c) is part of GLU.

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
   @file LU_SSE.c
   @brief LU decompositions
 */
#include "Mainfile.h"

// we have better determinant routines for small matrices
#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC ) && (NC>3)

#include <immintrin.h>
#include "LU.h"         // for alphabetising
#include "SSE2_OPS.h"   // SSE2_* macros

//  LU determinant computation with intrinsics
double complex
LU_det( const size_t N , const GLU_complex U[ N*N ] )
{
  GLU_complex a[ N*N ] ;
  memcpy( a , U , N*N * sizeof( GLU_complex ) ) ;
  return LU_det_overwrite( N , a ) ;
}

//  same as above overwrites U
double complex
LU_det_overwrite( const size_t N , GLU_complex U[ N*N ] )
{
  __m128d determinant = _mm_setr_pd( 1.0 , 0.0 ) ;
  register __m128d best , attempt , z1 , z2 ; 
  double complex s ;
  size_t i , j , l , piv , perms = 0 ;
  __m128d *a = (__m128d*)U ;

  for( i = 0 ; i < N-1 ; i++ ) {
    // z1 = a[i*(N+1)] , z2 = Im(z1),Re(z1)
    z1 = a[i*(N+1)] ;
    z2 = _mm_shuffle_pd( z1 , z1 , 1 ) ;
    best = _mm_add_pd( _mm_mul_pd( z1 , z1 ) ,
		       _mm_mul_pd( z2 , z2 ) ) ;
    piv = i ;
    // compare frob norms of other elements
    for( j = i+1 ; j < N ; j++ ) {
      z1 = a[i+j*N] ;
      z2 = _mm_shuffle_pd( z1 , z1 , 1 ) ;
      attempt = _mm_add_pd( _mm_mul_pd( z1 , z1 ) ,
			    _mm_mul_pd( z2 , z2 ) ) ;
      if( _mm_ucomilt_sd( best , attempt ) ) { 
	piv = j ; 
	best = attempt ; 
      }
    }
    if( piv != i ) {
      // physically swap rows
      __m128d *p1 = a + i*N , *p2 = a + piv * N ;
      for( l = 0 ; l < N ; l++ ) {
	z1 = *p1 ; *p1++ = *p2 ; *p2++ = z1 ;
      }
      perms++ ;
    }
    if( _mm_ucomile_sd( best , _mm_setzero_pd() ) ) {
      fprintf( stderr , "[DETERMINANT] LU  Singular Matrix!!!\n" ) ;
      return 0.0 ;
    }
    // perform gaussian elimination
    const __m128d dt = _mm_div_pd( SSE2_CONJ( a[ i*(N+1) ] ) ,
				   best ) ;

    for( j = N-1 ; j > i ; j-- ) { // go up in other column
      __m128d *pB = a + i + j*N ;
      register const __m128d fac1 = SSE2_MUL( *pB , dt ) ; pB++ ;
      // go along the row performing the subtraction, there is no point in
      // subtracting elements where we have determined the best pivot, just the
      // columns to the right of the pivot
      const __m128d *pA = a + i*(N+1) + 1 ;
      for( l = 0 ; l < N - i - 1 ; l++ ) {
	*pB = _mm_sub_pd( *pB , SSE2_MUL( fac1 , *pA ) ) , pB++ , pA++ ;
      }
    }
    determinant = SSE2_MUL( determinant , a[ i*(N+1) ] ) ;
  }
  determinant = SSE2_MUL( determinant , a[ N*N-1 ] ) ;
  _mm_store_pd( (void*)&s , determinant ) ;
  return perms&1 ? -s : s ;
}

#endif // NC < 3
