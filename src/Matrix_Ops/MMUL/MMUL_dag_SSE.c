/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (MMUL_dag_SSE.c) is part of GLU.

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
   @file MMUL_dag.c
   @brief matrix multiply \f$ a = b \times c^{\dagger}\f$ , where \f$ a,b,c \in NCxNC \f$ GLU_complex matrices 
 */

#include "Mainfile.h"

#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC )

#include <immintrin.h>
#include "SSE2_OPS.h"

// a = b * c^{\dagger}
#ifndef multab_dag
void 
multab_dag( GLU_complex a[ NCNC ] , 
	    const GLU_complex b[ NCNC ] , 
	    const GLU_complex c[ NCNC ] )
{
  __m128d *A = ( __m128d* )a ;
  const __m128d *B = ( const __m128d* )b ;
  const __m128d *C = ( const __m128d* )c ;
#if NC == 3
  // top row
  //a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) + b[2] * conj( c[2] ) ; 
  *( A + 0 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 0 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJ( *( B + 1 ) , *( C + 1 ) ) ,
				       SSE2_MUL_CONJ( *( B + 2 ) , *( C + 2 ) ) ) ) ;
  //a[1] = b[0] * conj( c[3] ) + b[1] * conj( c[4] ) + b[2] * conj( c[5] ) ; 
  *( A + 1 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 0 ) , *( C + 3 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJ( *( B + 1 ) , *( C + 4 ) ) ,
				       SSE2_MUL_CONJ( *( B + 2 ) , *( C + 5 ) ) ) ) ;
  //a[2] = b[0] * conj( c[6] ) + b[1] * conj( c[7] ) + b[2] * conj( c[8] ) ; 
  *( A + 2 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 0 ) , *( C + 6 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJ( *( B + 1 ) , *( C + 7 ) ) ,
				       SSE2_MUL_CONJ( *( B + 2 ) , *( C + 8 ) ) ) ) ;
  // middle row
  //a[3] = b[3] * conj( c[0] ) + b[4] * conj( c[1] ) + b[5] * conj( c[2] ) ; 
  *( A + 3 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 3 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJ( *( B + 4 ) , *( C + 1 ) ) ,
				       SSE2_MUL_CONJ( *( B + 5 ) , *( C + 2 ) ) ) ) ;

  //a[4] = b[3] * conj( c[3] ) + b[4] * conj( c[4] ) + b[5] * conj( c[5] ) ; 
  *( A + 4 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 3 ) , *( C + 3 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJ( *( B + 4 ) , *( C + 4 ) ) ,
				       SSE2_MUL_CONJ( *( B + 5 ) , *( C + 5 ) ) ) ) ;

  //a[5] = b[3] * conj( c[6] ) + b[4] * conj( c[7] ) + b[5] * conj( c[8] ) ;
  *( A + 5 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 3 ) , *( C + 6 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJ( *( B + 4 ) , *( C + 7 ) ) ,
				       SSE2_MUL_CONJ( *( B + 5 ) , *( C + 8 ) ) ) ) ;
  // bottom row
  //a[6] = b[6] * conj( c[0] ) + b[7] * conj( c[1] ) + b[8] * conj( c[2] ) ; 
  *( A + 6 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 6 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJ( *( B + 7 ) , *( C + 1 ) ) ,
				       SSE2_MUL_CONJ( *( B + 8 ) , *( C + 2 ) ) ) ) ;
  //a[7] = b[6] * conj( c[3] ) + b[7] * conj( c[4] ) + b[8] * conj( c[5] ) ; 
  *( A + 7 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 6 ) , *( C + 3 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJ( *( B + 7 ) , *( C + 4 ) ) ,
				       SSE2_MUL_CONJ( *( B + 8 ) , *( C + 5 ) ) ) ) ;
  //a[8] = b[6] * conj( c[6] ) + b[7] * conj( c[7] ) + b[8] * conj( c[8] ) ; 
  *( A + 8 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 6 ) , *( C + 6 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJ( *( B + 7 ) , *( C + 7 ) ) ,
				       SSE2_MUL_CONJ( *( B + 8 ) , *( C + 8 ) ) ) ) ;
#elif NC==2
  *( A + 0 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 0 ) , *( C + 0 ) ) ,
			   SSE2_MUL_CONJ( *( B + 1 ) , *( C + 1 ) ) ) ;
  *( A + 1 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 0 ) , *( C + 2 ) ) , 
			   SSE2_MUL_CONJ( *( B + 1 ) , *( C + 3 ) ) ) ;
  *( A + 2 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 2 ) , *( C + 0 ) ) ,
			   SSE2_MUL_CONJ( *( B + 3 ) , *( C + 1 ) ) ) ;
  *( A + 3 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 2 ) , *( C + 2 ) ) , 
			   SSE2_MUL_CONJ( *( B + 3 ) , *( C + 3 ) ) ) ;
#else // instead of inlining we have a function call
  size_t i , j , m ;

  // transpose C
  __m128d CT[ NCNC ] ; 
  for( i = 0 ; i < NC ; i++ ) {
    CT[ i*(NC+1) ] = C[ i*(NC+1)  ] ;
    for( j = i+1 ; j < NC ; j++ ) {
      CT[ j + i*NC ] = C[ i + j*NC ] ;
      CT[ i + j*NC ] = C[ j + i*NC ] ;
    }
    for( j = 0 ; j < NC ; j++ ) {
      *A = _mm_setzero_pd() ; A++ ;
    }
  }

  // loop in a weird order to help the cache
  // this is basically the MMUL_SSE routine 
  for( i = 0 ; i < NC ; i++ ) {
    C = (const __m128d*)CT ;
    for( m = 0 ; m < NC ; m++  ) {
      A = (__m128d*)( a + i*NC ) ;
      const __m128d mul = *B ; B++ ;
      #if (NC%2==1)
      *A = _mm_add_pd( *A , SSE2_MUL_CONJ( mul , *C ) ) ; A++ ; C++ ;
      for( j = 0 ; j < (NC-1)/NBLOCK ; j++ ) {
	M_REPEAT(NBLOCK,
		 *A = _mm_add_pd( *A , SSE2_MUL_CONJ( mul , *C ) ) ;\
		 A++ ; C++ ;)
      }
      #else
      for( j = 0 ; j < NC/NBLOCK ; j++ ) {
	M_REPEAT(NBLOCK,
		 *A = _mm_add_pd( *A , SSE2_MUL_CONJ( mul , *C ) ) ;\
		 A++ ; C++ ;)
      }
      #endif
    }
  }
#endif
  return ;
}
#endif //!defined multab_dag

#endif // end of IMMINTRIN_H
