/*
    Copyright 2013 Renwick James Hudspith

    This file (MMULdagdag.c) is part of GLU.

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
   @file MMULdagdag.c
   @brief computes \f$ a = b^{\dagger} c^{\dagger} \f$ matrix multiplication
 */

#include "Mainfile.h"

#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC )

#include <immintrin.h>
#include "SSE2_OPS.h"

#ifndef multab_dagdag
// a = b^{\dagger} * c^{\dagger}
void 
multab_dagdag( GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC ] , 
	       const GLU_complex c[ NCNC ] )
{
  // recast to alignment
  __m128d *A = (__m128d*)a ;
  const __m128d *B = (const __m128d*)b ;
  const __m128d *C = (const __m128d*)c ;
#if NC==3
  // a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] ) + conj( b[6] ) * conj( c[2] ) ;
  *( A + 0 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 3 ) , *( C + 1 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 6 ) , *( C + 2 ) ) ) ) ;
  //a[1] = conj( b[0] ) * conj( c[3] ) + conj( b[3] ) * conj( c[4] ) + conj( b[6] ) * conj( c[5] ) ; 
  *( A + 1 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 3 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 3 ) , *( C + 4 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 6 ) , *( C + 5 ) ) ) ) ;
  //a[2] = conj( b[0] ) * conj( c[6] ) + conj( b[3] ) * conj( c[7] ) + conj( b[6] ) * conj( c[8] ) ; 
  *( A + 2 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 6 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 3 ) , *( C + 7 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 6 ) , *( C + 8 ) ) ) ) ;
  // middle row //
  //a[3] = conj( b[1] ) * conj( c[0] ) + conj( b[4] ) * conj( c[1] ) + conj( b[7] ) * conj( c[2] ) ; 
  *( A + 3 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 1 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 4 ) , *( C + 1 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 7 ) , *( C + 2 ) ) ) ) ;
  //a[4] = conj( b[1] ) * conj( c[3] ) + conj( b[4] ) * conj( c[4] ) + conj( b[7] ) * conj( c[5] ) ; 
  *( A + 4 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 1 ) , *( C + 3 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 4 ) , *( C + 4 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 7 ) , *( C + 5 ) ) ) ) ;
  //a[5] = conj( b[1] ) * conj( c[6] ) + conj( b[4] ) * conj( c[7] ) + conj( b[7] ) * conj( c[8] ) ; 
  *( A + 5 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 1 ) , *( C + 6 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 4 ) , *( C + 7 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 7 ) , *( C + 8 ) ) ) ) ;
  //a[6] = conj( b[2] ) * conj( c[0] ) + conj( b[5] ) * conj( c[1] ) + conj( b[8] ) * conj( c[2] ) ;
  *( A + 6 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 2 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 5 ) , *( C + 1 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 8 ) , *( C + 2 ) ) ) ) ;
  //a[7] = conj( b[2] ) * conj( c[3] ) + conj( b[5] ) * conj( c[4] ) + conj( b[8] ) * conj( c[5] ) ;
  *( A + 7 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 2 ) , *( C + 3 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 5 ) , *( C + 4 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 8 ) , *( C + 5 ) ) ) ) ;
  //a[8] = conj( b[2] ) * conj( c[6] ) + conj( b[5] ) * conj( c[7] ) + conj( b[8] ) * conj( c[8] ) ; 
  *( A + 8 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 2 ) , *( C + 6 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 5 ) , *( C + 7 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 8 ) , *( C + 8 ) ) ) ) ;
#elif NC==2
  //a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[2] ) * conj( c[1] ) ;
  *( A + 0 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 0 ) ) ,
			   SSE2_MUL_CONJCONJ( *( B + 2 ) , *( C + 1 ) ) ) ;
  //a[1] = conj( b[0] ) * conj( c[2] ) + conj( b[2] ) * conj( c[3] ) ;
  *( A + 1 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 2 ) ) ,
			   SSE2_MUL_CONJCONJ( *( B + 2 ) , *( C + 3 ) ) ) ;
  //a[2] = conj( b[1] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] ) ; 
  *( A + 2 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 1 ) , *( C + 0 ) ) ,
			   SSE2_MUL_CONJCONJ( *( B + 3 ) , *( C + 1 ) ) ) ;
  //a[3] = conj( b[1] ) * conj( c[2] ) + conj( b[3] ) * conj( c[3] ) ; 
  *( A + 3 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 1 ) , *( C + 2 ) ) ,
			   SSE2_MUL_CONJCONJ( *( B + 3 ) , *( C + 3 ) ) ) ;
#else
  int i , j , m ;
  register GLU_complex sumr = 0.0 , sumi = 0.0 ;
  register GLU_real REB , IMB , REC , IMC ;
  const GLU_complex *pC , *pB ;
  for( i = 0 ; i < NC ; i++ ) {
    pC = c ;
    for( j = 0 ; j < NC ; j++ ) {
      pB = b ;
      sumr = sumi = 0.0 ;
      for( m = 0 ; m < NC ; m++ ) {
	REB = creal( pB[i] ) ; IMB = cimag( pB[i] ) ;
	REC = creal( pC[m] ) ; IMC = cimag( pC[m] ) ;
	sumr += REB * REC - IMB * IMC ;
	sumi += REB * IMC + IMB * REC ;
	pB += NC ;
      }
      a[ j + NC*i ] = sumr + I * sumi ;
      pC += NC ;
    }
  }
#endif
  return ;
}
#endif

#ifndef multab_dagdag_suNC

void //__attribute__((hot))
multab_dagdag_suNC( GLU_complex a[ NCNC ] , 
		    const GLU_complex b[ NCNC ] , 
		    const GLU_complex c[ NCNC ] )
{
  // recast to alignment
  __m128d *A = (__m128d*)a ;
  const __m128d *B = (const __m128d*)b ;
  const __m128d *C = (const __m128d*)c ;
#if NC == 3
  // top row //
  // a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] ) + conj( b[6] ) * conj( c[2] ) ;
  *( A + 0 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 3 ) , *( C + 1 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 6 ) , *( C + 2 ) ) ) ) ;
  //a[1] = conj( b[0] ) * conj( c[3] ) + conj( b[3] ) * conj( c[4] ) + conj( b[6] ) * conj( c[5] ) ; 
  *( A + 1 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 3 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 3 ) , *( C + 4 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 6 ) , *( C + 5 ) ) ) ) ;
  //a[2] = conj( b[0] ) * conj( c[6] ) + conj( b[3] ) * conj( c[7] ) + conj( b[6] ) * conj( c[8] ) ; 
  *( A + 2 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 6 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 3 ) , *( C + 7 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 6 ) , *( C + 8 ) ) ) ) ;
  // middle row //
  //a[3] = conj( b[1] ) * conj( c[0] ) + conj( b[4] ) * conj( c[1] ) + conj( b[7] ) * conj( c[2] ) ; 
  *( A + 3 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 1 ) , *( C + 0 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 4 ) , *( C + 1 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 7 ) , *( C + 2 ) ) ) ) ;
  //a[4] = conj( b[1] ) * conj( c[3] ) + conj( b[4] ) * conj( c[4] ) + conj( b[7] ) * conj( c[5] ) ; 
  *( A + 4 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 1 ) , *( C + 3 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 4 ) , *( C + 4 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 7 ) , *( C + 5 ) ) ) ) ;
  //a[5] = conj( b[1] ) * conj( c[6] ) + conj( b[4] ) * conj( c[7] ) + conj( b[7] ) * conj( c[8] ) ; 
  *( A + 5 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 1 ) , *( C + 6 ) ) ,
			   _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 4 ) , *( C + 7 ) ) ,
				       SSE2_MUL_CONJCONJ( *( B + 7 ) , *( C + 8 ) ) ) ) ;
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
  //a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[2] ) * conj( c[1] ) ;
  *( A + 0 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 0 ) ) ,
			   SSE2_MUL_CONJCONJ( *( B + 2 ) , *( C + 1 ) ) ) ;
  //a[1] = conj( b[0] ) * conj( c[2] ) + conj( b[2] ) * conj( c[3] ) ;
  *( A + 1 ) = _mm_add_pd( SSE2_MUL_CONJCONJ( *( B + 0 ) , *( C + 2 ) ) ,
			   SSE2_MUL_CONJCONJ( *( B + 2 ) , *( C + 3 ) ) ) ;
  // a[2] = -conj( a[1] )
  *( A + 2 ) = SSE2_FLIP( SSE2_CONJ( *( A + 1 ) ) ) ; 
  // a[3] =  conj( a[0] )
  *( A + 3 ) = SSE2_CONJ( *( A + 0 ) ) ;
#else
  return multab_dagdag( a , b , c ) ;
#endif
  return ;
}
#endif

#endif // immintrin.h