/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (MMUL.c) is part of GLU.

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
   @file MMUL_SUNC_SSE.c
   @brief SU(NC) matrix-matrix multiply with SSE intrinsics
 */

#include "Mainfile.h"

#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC )  

#ifndef multab_suNC

#include <immintrin.h>
#include "SSE2_OPS.h"

void
multab_suNC( GLU_complex a[ NCNC ] , 
	     const GLU_complex b[ NCNC ] , 
	     const GLU_complex c[ NCNC ] )
{
 // recast to alignment
  __m128d *A = (__m128d*)a ;
  const __m128d *B = (const __m128d*)b ;
  const __m128d *C = (const __m128d*)c ;
#if NC == 3
  // top row
  //a[0] = b[0] * c[0] + b[1] * c[3] + b[2] * c[6] ; 
  *( A + 0 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( C + 3 ) ) ,
				       SSE2_MUL( *( B + 2 ) , *( C + 6 ) ) ) ) ;
  // a[1] = b[0] * c[1] + b[1] * c[4] + b[2] * c[7] ; 
  *( A + 1 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 1 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( C + 4 ) ) ,
				       SSE2_MUL( *( B + 2 ) , *( C + 7 ) ) ) ) ;
  // a[2] = b[0] * c[2] + b[1] * c[5] + b[2] * c[8] ;
  *( A + 2 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 2 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 1 ) , *( C + 5 ) ) ,
				       SSE2_MUL( *( B + 2 ) , *( C + 8 ) ) ) ) ;
  // middle row //
  // a[3] = b[3] * c[0] + b[4] * c[3] + b[5] * c[6] ; 
  *( A + 3 ) = _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( C + 3 ) ) ,
				       SSE2_MUL( *( B + 5 ) , *( C + 6 ) ) ) ) ;
  // a[4] = b[3] * c[1] + b[4] * c[4] + b[5] * c[7] ; 
  *( A + 4 ) = _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( C + 1 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( C + 4 ) ) ,
				       SSE2_MUL( *( B + 5 ) , *( C + 7 ) ) ) ) ;
  // a[5] = b[3] * c[2] + b[4] * c[5] + b[5] * c[8] ; 
  *( A + 5 ) = _mm_add_pd( SSE2_MUL( *( B + 3 ) , *( C + 2 ) ) ,
			   _mm_add_pd( SSE2_MUL( *( B + 4 ) , *( C + 5 ) ) ,
				       SSE2_MUL( *( B + 5 ) , *( C + 8 ) ) ) ) ;
  // bottom row // as a completion of the top two
  // a[6] = conj( a[1] * a[5] - a[2] * a[4] ) ; 
  *( A + 6 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( A + 1 ) , *( A + 5 ) ) ,
				      SSE2_MUL( *( A + 2 ) , *( A + 4 ) ) ) ) ;
  // a[7] = conj( a[2] * a[3] - a[0] * a[5] ) ; 
  *( A + 7 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( A + 2 ) , *( A + 3 ) ) ,
				      SSE2_MUL( *( A + 0 ) , *( A + 5 ) ) ) ) ;
  // a[8] = conj( a[0] * a[4] - a[1] * a[3] ) ; 
  *( A + 8 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( A + 0 ) , *( A + 4 ) ) ,
				      SSE2_MUL( *( A + 1 ) , *( A + 3 ) ) ) ) ;
#elif NC == 2
  *( A + 0 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 0 ) ) ,
			   SSE2_MUL( *( B + 1 ) , *( C + 2 ) ) ) ;
  *( A + 1 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , *( C + 1 ) ) ,
			   SSE2_MUL( *( B + 1 ) , *( C + 3 ) ) ) ;
  // a[2] = -conj( a[1] )
  *( A + 2 ) = SSE_FLIP( SSE2_CONJ( *( A + 1 ) ) ) ; 
  // a[3] =  conj( a[0] )
  *( A + 3 ) = SSE2_CONJ( *( A + 0 ) ) ;
#else
  return multab_suNC( a , b , c ) ;
#endif
  return ;
}

#endif // ifndef multab_suNC

#endif // end of HAVE_IMMINTRIN_H

