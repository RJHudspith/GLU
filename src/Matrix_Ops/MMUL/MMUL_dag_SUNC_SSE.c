/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (MMUL_dag_SUNC_SSE.c) is part of GLU.

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
   @file MMUL_dag_SUNC_SSE.c
   @brief matrix multiply \f$ a = b \times c^{\dagger}\f$ , where \f$ a,b,c \in NCxNC \f$ GLU_complex matrices 
 */
#include "Mainfile.h"

#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC )

#include <immintrin.h>
#include "SSE2_OPS.h"

#ifndef multab_dag_suNC

void //__attribute__((hot))
multab_dag_suNC( GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] , 
		 const GLU_complex c[ NCNC ] )
{
  __m128d *A = ( __m128d* )a ;
  const __m128d *B = ( const __m128d* )b ;
  const __m128d *C = ( const __m128d* )c ;
#if NC == 3
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
  //a[0] = b[0]*conj( c[0] ) + b[1]*conj( c[1] )
  *( A + 0 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 0 ) , *( C + 0 ) ) ,
			   SSE2_MUL_CONJ( *( B + 1 ) , *( C + 1 ) ) ) ;
  //a[1] = b[0] * conj( c[2] ) + b[1] * conj( c[3] ) 
  *( A + 1 ) = _mm_add_pd( SSE2_MUL_CONJ( *( B + 0 ) , *( C + 2 ) ) , 
			   SSE2_MUL_CONJ( *( B + 1 ) , *( C + 3 ) ) ) ;
  // a[2] = -conj( a[1] )
  *( A + 2 ) = SSE_FLIP( SSE2_CONJ( *( A + 1 ) ) ) ; 
  // a[3] =  conj( a[0] )
  *( A + 3 ) = SSE2_CONJ( *( A + 0 ) ) ;
#else
  return multab_dag( a , b , c ) ;
#endif
  return ;
}
#endif

#endif // end of IMMINTRIN_H
