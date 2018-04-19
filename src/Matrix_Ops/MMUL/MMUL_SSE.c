/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (MMUL_SSE.c) is part of GLU.

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
   @file MMUL_SSE.c
   @brief matrix-matrix multiply with SSE intrinsics
 */
#include "Mainfile.h"

#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC )  

#include <immintrin.h>
#include "SSE2_OPS.h"

#ifndef multab_atomic_left
// simple matrix multiplication ( left multiply ) a = b * a 
void 
multab_atomic_left( GLU_complex a[ NCNC ] , 
		    const GLU_complex b[ NCNC ] )
{
  __m128d *A = ( __m128d* )a ;
  const __m128d *B = ( const __m128d* )b ;
#if NC == 3
  register __m128d C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 0 ) ) , 
				    _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( A + 3 ) ) ,
						SSE2_MUL( *( B + 2 ) , *( A + 6 ) ) ) ) ;
  register __m128d C1 = _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( A + 0 ) ) , 
				    _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( A + 3 ) ) ,
						SSE2_MUL( *( B + 5 ) , *( A + 6 ) ) ) ) ;
  register __m128d C2 = _mm_add_pd( SSE2_MUL( *( B + 6 ) , *( A + 0 ) ) , 
				    _mm_add_pd( SSE2_MUL( *( B + 7 ) , *( A + 3 ) ) ,
						SSE2_MUL( *( B + 8 ) , *( A + 6 ) ) ) ) ;
  A[0] = C0 ; A[3] = C1 ; A[6] = C2 ;
  // middle
  C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 1 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( A + 4 ) ) ,
			       SSE2_MUL( *( B + 2 ) , *( A + 7 ) ) ) ) ;
  C1 = _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( A + 1 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( A + 4 ) ) ,
			       SSE2_MUL( *( B + 5 ) , *( A + 7 ) ) ) ) ;
  C2 = _mm_add_pd( SSE2_MUL( *( B + 6 ) , *( A + 1 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 7 ) , *( A + 4 ) ) ,
			       SSE2_MUL( *( B + 8 ) , *( A + 7 ) ) ) ) ;
  A[1] = C0 ; A[4] = C1 ; A[7] = C2 ;
  // end
  C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 2 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( A + 5 ) ) ,
			       SSE2_MUL( *( B + 2 ) , *( A + 8 ) ) ) ) ;
  C1 = _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( A + 2 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( A + 5 ) ) ,
			       SSE2_MUL( *( B + 5 ) , *( A + 8 ) ) ) ) ;
  C2 = _mm_add_pd( SSE2_MUL( *( B + 6 ) , *( A + 2 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 7 ) , *( A + 5 ) ) ,
			       SSE2_MUL( *( B + 8 ) , *( A + 8 ) ) ) ) ;
  A[2] = C0 ; A[5] = C1 ; A[8] = C2 ;
#elif NC == 2
  __m128d C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 0 ) ) ,
			   SSE2_MUL( *( B + 1 ) , *( A + 2 ) ) ) ;
  __m128d C1 = _mm_add_pd( SSE2_MUL( *( B + 2 ) , *( A + 0 ) ) ,
			   SSE2_MUL( *( B + 3 ) , *( A + 2 ) ) ) ;
  A[0] = C0 ; A[2] = C1 ;
  C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 1 ) ) ,
		   SSE2_MUL( *( B + 1 ) , *( A + 3 ) ) ) ;
  C1 = _mm_add_pd( SSE2_MUL( *( B + 2 ) , *( A + 1 ) ) ,
		   SSE2_MUL( *( B + 3 ) , *( A + 3 ) ) ) ;
  A[1] = C0 ; A[3] = C1 ;					
#else
  // slow and stupid loopy version, memory access pattern is odd
  size_t i , j , m ;
  __m128d C[ NC ] ; // temporary storage up in here
  register __m128d sum ;
  for( i = 0 ; i < NC ; i++ ) { // loop cols
    B = (const __m128d*)b ;
    for( j = 0 ; j < NC ; j ++ ) { // loop rows
      sum = _mm_setzero_pd( ) ;
      for( m = 0 ; m < NC ; m ++  ) { // loop elements in row or column
	sum = _mm_add_pd( sum , SSE2_MUL( *( B ) , *( A + i + m * NC ) ) ) ;
	B++ ;
      }	
      *( C + j ) = sum ;
    }
    // copy back over to a ...
    for( m = 0 ; m < NC ; m++ ) {
      *( A + i + m*NC ) = *( C + m ) ;
    }
  }
#endif
  return ;
}
#endif

#ifndef multab_atomic_right
// simple matrix multiplication a = a * b
void 
multab_atomic_right( GLU_complex a[ NCNC ] , 
		     const GLU_complex b[ NCNC ] )
{
  __m128d *A = ( __m128d* )a ;
  const __m128d *B = ( const __m128d* )b ;
#if NC==3
  register __m128d C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 0 ) ) , 
				    _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( A + 1 ) ) ,
						SSE2_MUL( *( B + 6 ) , *( A + 2 ) ) ) ) ;
  register __m128d C1 = _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( A + 0 ) ) , 
				    _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( A + 1 ) ) ,
						SSE2_MUL( *( B + 7 ) , *( A + 2 ) ) ) ) ;
  register __m128d C2 = _mm_add_pd( SSE2_MUL( *( B + 2 ) , *( A + 0 ) ) , 
				    _mm_add_pd( SSE2_MUL( *( B + 5 ) , *( A + 1 ) ) ,
						SSE2_MUL( *( B + 8 ) , *( A + 2 ) ) ) ) ;
  A[0] = C0 ; A[1] = C1 ; A[2] = C2 ;
  // middle
  C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 3 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( A + 4 ) ) ,
			       SSE2_MUL( *( B + 6 ) , *( A + 5 ) ) ) ) ;
  C1 = _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( A + 3 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( A + 4 ) ) ,
			       SSE2_MUL( *( B + 7 ) , *( A + 5 ) ) ) ) ;
  C2 = _mm_add_pd( SSE2_MUL( *( B + 2 ) , *( A + 3 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 5 ) , *( A + 4 ) ) ,
			       SSE2_MUL( *( B + 8 ) , *( A + 5 ) ) ) ) ;
  A[3] = C0 ; A[4] = C1 ; A[5] = C2 ;
  // bottom
  C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 6 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( A + 7 ) ) ,
			       SSE2_MUL( *( B + 6 ) , *( A + 8 ) ) ) ) ;
  C1 = _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( A + 6 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( A + 7 ) ) ,
			       SSE2_MUL( *( B + 7 ) , *( A + 8 ) ) ) ) ;
  C2 = _mm_add_pd( SSE2_MUL( *( B + 2 ) , *( A + 6 ) ) , 
		   _mm_add_pd( SSE2_MUL( *( B + 5 ) , *( A + 7 ) ) ,
			       SSE2_MUL( *( B + 8 ) , *( A + 8 ) ) ) ) ;
  A[6] = C0 ; A[7] = C1 ; A[8] = C2 ;				
#elif NC==2
  __m128d C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 0 ) ) ,
			   SSE2_MUL( *( B + 2 ) , *( A + 1 ) ) ) ;
  __m128d C1 = _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( A + 0 ) ) ,
			   SSE2_MUL( *( B + 3 ) , *( A + 1 ) ) ) ;
  A[0] = C0 ; A[1] = C1 ;
  C0 = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( A + 2 ) ) ,
		   SSE2_MUL( *( B + 2 ) , *( A + 3 ) ) ) ;
  C1 = _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( A + 2 ) ) ,
		   SSE2_MUL( *( B + 3 ) , *( A + 3 ) ) ) ;
  A[2] = C0 ; A[3] = C1 ;				
#else
  // slow and stupid loopy version
  size_t i , j , m ;
  __m128d C[ NC ] ;
  register __m128d sum ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j ++ ) {
      sum = _mm_setzero_pd( ) ;
      for( m = 0 ; m < NC ; m ++  ) {
        sum = _mm_add_pd( sum , 
			  SSE2_MUL( *( A +  m + i*NC ) , *( B + j + m*NC ) ) ) ;
      }
      *( C + j ) = sum ;
    }
    // copy back over to a ...
    for( j = 0 ; j < NC ; j ++ ) {
      *( A + j + i*NC ) = *( C + j ) ;
    }
  }
#endif
  return ;
}
#endif

#ifndef multab
// simple matrix multiplication
void 
multab( GLU_complex a[ NCNC ] , 
	const GLU_complex b[ NCNC ] , 
	const GLU_complex c[ NCNC ] )
{
  __m128d *A = ( __m128d* )a ;
  const __m128d *B = ( __m128d* )b ;
  const __m128d *C = ( __m128d* )c ;
#if NC==3
  *( A + 0 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( C + 3 ) ) ,
				       SSE2_MUL( *( B + 2 ) , *( C + 6 ) ) ) ) ;
  *( A + 1 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 1 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( C + 4 ) ) ,
				       SSE2_MUL( *( B + 2 ) , *( C + 7 ) ) ) ) ;
  *( A + 2 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 2 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( C + 5 ) ) ,
				       SSE2_MUL( *( B + 2 ) , *( C + 8 ) ) ) ) ;
  *( A + 3 ) = _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( C + 3 ) ) ,
				       SSE2_MUL( *( B + 5 ) , *( C + 6 ) ) ) ) ;
  *( A + 4 ) = _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( C + 1 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( C + 4 ) ) ,
				       SSE2_MUL( *( B + 5 ) , *( C + 7 ) ) ) ) ;
  *( A + 5 ) = _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( C + 2 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( C + 5 ) ) ,
				       SSE2_MUL( *( B + 5 ) , *( C + 8 ) ) ) ) ;
  *( A + 6 ) = _mm_add_pd( SSE2_MUL( *( B + 6 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 7 ) , *( C + 3 ) ) ,
				       SSE2_MUL( *( B + 8 ) , *( C + 6 ) ) ) ) ;
  *( A + 7 ) = _mm_add_pd( SSE2_MUL( *( B + 6 ) , *( C + 1 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 7 ) , *( C + 4 ) ) ,
				       SSE2_MUL( *( B + 8 ) , *( C + 7 ) ) ) ) ;
  *( A + 8 ) = _mm_add_pd( SSE2_MUL( *( B + 6 ) , *( C + 2 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 7 ) , *( C + 5 ) ) ,
				       SSE2_MUL( *( B + 8 ) , *( C + 8 ) ) ) ) ;
#elif NC==2
  *( A + 0 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 0 ) ) ,
			   SSE2_MUL( *( B + 1 ) , *( C + 2 ) ) ) ;
  *( A + 1 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 1 ) ) ,
			   SSE2_MUL( *( B + 1 ) , *( C + 3 ) ) ) ;
  *( A + 2 ) = _mm_add_pd( SSE2_MUL( *( B + 2 ) , *( C + 0 ) ) ,
			   SSE2_MUL( *( B + 3 ) , *( C + 2 ) ) ) ;
  *( A + 3 ) = _mm_add_pd( SSE2_MUL( *( B + 2 ) , *( C + 1 ) ) ,
			   SSE2_MUL( *( B + 3 ) , *( C + 3 ) ) ) ;
#else
  // standard 3 loop version
  size_t i , j , m ;
  // init the result matrix "a" to zero
  for( i = 0 ; i < NCNC ; i++ ) {
    *A = _mm_setzero_pd( ) , A++ ;
  }
  // loop in a weird order to help the cache
  for( i = 0 ; i < NC ; i++ ) {
    C = (__m128d*)( c ) ;
    for( m = 0 ; m < NC ; m++  ) {
      A = (__m128d*)( a + i*NC ) ;
      const __m128d mul = *B ; B++ ;
      #if (NC%2==1)
      *A = _mm_add_pd( *A , SSE2_MUL( mul , *C ) ) ; A++ ; C++ ;
      for( j = 0 ; j < (NC-1)/NBLOCK ; j++ ) {
	M_REPEAT(NBLOCK,
		 *A = _mm_add_pd( *A , SSE2_MUL( mul , *C ) ) ;\
		 A++ ; C++ ;)
      }
      #else
      for( j = 0 ; j < NC/NBLOCK ; j++ ) {
	M_REPEAT(NBLOCK,
		 *A = _mm_add_pd( *A , SSE2_MUL( mul , *C ) ) ;\
		 A++ ; C++ ;)
      }
      #endif
    }
  }
#endif
  return ;
}
#endif

#endif // end of HAVE_IMMINTRIN_H

