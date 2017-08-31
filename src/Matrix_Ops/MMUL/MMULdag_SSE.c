/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (MMULdag.c) is part of GLU.

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
   @file MMULdag.c
   @brief computes \f$ a = b^{\dagger}\times c \f$
 */

#include "Mainfile.h"

#if (defined HAVE_IMMINTRIN_H) && !( defined SINGLE_PREC )

#include <immintrin.h>
#include "SSE2_OPS.h"

#ifndef multabdag
// 3x3 mult a=b^{\dagger}.c//
void 
multabdag( GLU_complex a[ NCNC ] , 
	   const GLU_complex b[ NCNC ] , 
	   const GLU_complex c[ NCNC ] )
{
  // recast to alignment
  __m128d *A = (__m128d*)a ;
  const __m128d *B = (const __m128d*)b ;
  const __m128d *C = (const __m128d*)c ;
#if NC==3
  // a[0] = conj( b[0] ) * c[0] + conj( b[3] ) * c[3] + conj( b[6] ) * c[6] ;
  *( A + 0 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 0 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( *( B + 3 ) , *( C + 3 ) ) ,
				       SSE2_MULCONJ( *( B + 6 ) , *( C + 6 ) ) ) ) ;
  //a[1] = conj( b[0] ) * c[1] + conj( b[3] ) * c[4] + conj( b[6] ) * c[7] ;
  *( A + 1 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 0 ) , *( C + 1 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( *( B + 3 ) , *( C + 4 ) ) ,
				       SSE2_MULCONJ( *( B + 6 ) , *( C + 7 ) ) ) ) ;
  //a[2] = conj( b[0] ) * c[2] + conj( b[3] ) * c[5] + conj( b[6] ) * c[8] ;
  *( A + 2 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 0 ) , *( C + 2 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( *( B + 3 ) , *( C + 5 ) ) ,
				       SSE2_MULCONJ( *( B + 6 ) , *( C + 8 ) ) ) ) ;
  // middle row
  //a[3] = conj( b[1] ) * c[0] + conj( b[4] ) * c[3] + conj( b[7] ) * c[6] ;
  *( A + 3 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 1 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( *( B + 4 ) , *( C + 3 ) ) ,
				       SSE2_MULCONJ( *( B + 7 ) , *( C + 6 ) ) ) ) ;
  //a[4] = conj( b[1] ) * c[1] + conj( b[4] ) * c[4] + conj( b[7] ) * c[7] ;
  *( A + 4 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 1 ) , *( C + 1 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( *( B + 4 ) , *( C + 4 ) ) ,
				       SSE2_MULCONJ( *( B + 7 ) , *( C + 7 ) ) ) ) ;
  //a[5] = conj( b[1] ) * c[2] + conj( b[4] ) * c[5] + conj( b[7] ) * c[8] ;
  *( A + 5 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 1 ) , *( C + 2 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( *( B + 4 ) , *( C + 5 ) ) ,
				       SSE2_MULCONJ( *( B + 7 ) , *( C + 8 ) ) ) ) ;
  //a[6] = conj( b[2] ) * c[0] + conj( b[5] ) * c[3] + conj( b[8] ) * c[6] ; 
  *( A + 6 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 2 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( *( B + 5 ) , *( C + 3 ) ) ,
				       SSE2_MULCONJ( *( B + 8 ) , *( C + 6 ) ) ) ) ;
  //a[7] = conj( b[2] ) * c[1] + conj( b[5] ) * c[4] + conj( b[8] ) * c[7] ; 
  *( A + 7 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 2 ) , *( C + 1 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( *( B + 5 ) , *( C + 4 ) ) ,
				       SSE2_MULCONJ( *( B + 8 ) , *( C + 7 ) ) ) ) ;
  //a[8] = conj( b[2] ) * c[2] + conj( b[5] ) * c[5] + conj( b[8] ) * c[8] ;
  *( A + 8 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 2 ) , *( C + 2 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( *( B + 5 ) , *( C + 5 ) ) ,
				       SSE2_MULCONJ( *( B + 8 ) , *( C + 8 ) ) ) ) ;
#elif NC==2
  //a[0] = conj( b[0] ) * c[0] + conj( b[2] ) * c[2] ;
  *( A + 0 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 0 ) , *( C + 0 ) ) ,
			   SSE2_MULCONJ( *( B + 2 ) , *( C + 2 ) ) ) ;
  //a[1] = conj( b[0] ) * c[1] + conj( b[2] ) * c[3] ;
  *( A + 1 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 0 ) , *( C + 1 ) ) ,
			   SSE2_MULCONJ( *( B + 2 ) , *( C + 3 ) ) ) ;
  //a[2] = conj( b[1] ) * c[0] + conj( b[3] ) * c[2] ;
  *( A + 2 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 1 ) , *( C + 0 ) ) ,
			   SSE2_MULCONJ( *( B + 3 ) , *( C + 2 ) ) ) ;
  //a[3] = conj( b[1] ) * c[1] + conj( b[3] ) * c[3] ;
  *( A + 3 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 1 ) , *( C + 1 ) ) ,
			   SSE2_MULCONJ( *( B + 3 ) , *( C + 3 ) ) ) ;
#else
  size_t i , j , m ;
  for( i = 0 ; i < NCNC ; i++ ) {
    *A = _mm_setzero_pd() ; A++ ;
  }
  // loop designed for good cache througput
  for( m = 0 ; m < NC ; m++ ) {
    A = (__m128d*)a ;
    for( i = 0 ; i < NC ; i++ ) {
      C = (const __m128d*)( c + NC*m ) ;
      register const __m128d mul = *B ; B++ ;
      #if (NC%2==1)
      *A = _mm_add_pd( *A , SSE2_MULCONJ( mul , *C ) ) ; A++ ; C++ ;
      for( j = 0 ; j < (NC-1)/NBLOCK ; j++ ) {
	M_REPEAT(NBLOCK,
		 *A = _mm_add_pd( *A , SSE2_MULCONJ( mul , *C ) ) ;	\
		 A++ ; C++ ;)
      }
      #else
      for( j = 0 ; j < NC/NBLOCK ; j++ ) {
	M_REPEAT(NBLOCK,
		 *A = _mm_add_pd( *A , SSE2_MULCONJ( mul , *C ) ) ;	\
		 A++ ; C++ ;)
	  }
      #endif
    }
  }
#endif
  return ;
}
#endif // !defined multabdag

#endif
