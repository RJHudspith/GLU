/*
    Copyright 2013 Renwick James Hudspith

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
   @file MMUL.c
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
#if NC == 3
  GLU_complex C0 = b[0] * a[0] + b[1] * a[3] + b[2] * a[6] ;	\
  GLU_complex C1 = b[3] * a[0] + b[4] * a[3] + b[5] * a[6] ;	\
  GLU_complex C2 = b[6] * a[0] + b[7] * a[3] + b[8] * a[6] ;	\
  a[0] = C0 ; a[3] = C1 ; a[6] = C2 ;				\
  C0 = b[0] * a[1] + b[1] * a[4] + b[2] * a[7] ;		\
  C1 = b[3] * a[1] + b[4] * a[4] + b[5] * a[7] ;		\
  C2 = b[6] * a[1] + b[7] * a[4] + b[8] * a[7] ;		\
  a[1] = C0 ; a[4] = C1 ; a[7] = C2 ;				\
  C0 = b[0] * a[2] + b[1] * a[5] + b[2] * a[8] ;		\
  C1 = b[3] * a[2] + b[4] * a[5] + b[5] * a[8] ;		\
  C2 = b[6] * a[2] + b[7] * a[5] + b[8] * a[8] ;		\
  a[2] = C0 ; a[5] = C1 ; a[8] = C2 ;
#elif NC == 2
  GLU_complex C0 = b[0] * a[0] + b[1] * a[2] ;			\
  GLU_complex C1 = b[2] * a[0] + b[3] * a[2] ;			\
  a[0] = C0 ; a[2] = C1 ;					\
  C0 = b[0] * a[1] + b[1] * a[3] ;				\
  C1 = b[2] * a[1] + b[3] * a[3] ;				\
  a[1] = C0 ; a[3] = C1 ;					
#else
  // slow and stupid loopy version, memory access pattern is odd
  int i , j , m ;
  GLU_complex R[ NC ] ;
  register GLU_complex sum ;
  for( i = 0 ; i < NC ; i++ ) { // loop cols
    for( j = 0 ; j < NC ; j ++ ) { // loop rows
      sum = 0.0 ;
      for( m = 0 ; m < NC ; m ++  ) { // loop elements in row or column
        sum += b[ m + j*NC ] * a[ i + m*NC ] ;
      }	
      R[j] = sum ;
    }
    // copy back over to a ...
    for( m = 0 ; m < NC ; m++ ) {
      a[ i + m*NC ] = R[m] ;
    }
  }
#endif
  return ;
}
#endif

#ifndef multab_atomic
// simple matrix multiplication a = a * b
void 
multab_atomic_right( GLU_complex a[ NCNC ] , 
		     const GLU_complex b[ NCNC ] )
{
#if NC==3
  GLU_complex R0 = a[0] * b[0] + a[1] * b[3] + a[2] * b[6] ;	\
  GLU_complex R1 = a[0] * b[1] + a[1] * b[4] + a[2] * b[7] ;	\
  GLU_complex R2 = a[0] * b[2] + a[1] * b[5] + a[2] * b[8] ;	\
  a[0] = R0 ; a[1] = R1 ; a[2] = R2 ;				\
  R0 = a[3] * b[0] + a[4] * b[3] + a[5] * b[6] ;		\
  R1 = a[3] * b[1] + a[4] * b[4] + a[5] * b[7] ;		\
  R2 = a[3] * b[2] + a[4] * b[5] + a[5] * b[8] ;		\
  a[3] = R0 ; a[4] = R1 ; a[5] = R2 ;				\
  R0 = a[6] * b[0] + a[7] * b[3] + a[8] * b[6] ;		\
  R1 = a[6] * b[1] + a[7] * b[4] + a[8] * b[7] ;		\
  R2 = a[6] * b[2] + a[7] * b[5] + a[8] * b[8] ;		\
  a[6] = R0 ; a[7] = R1 ; a[8] = R2 ;				
#elif NC==2
  GLU_complex R0 = a[0] * b[0] + a[1] * b[2] ;		\
  GLU_complex R1 = a[0] * b[1] + a[1] * b[3] ;		\
  a[0] = R0 ; a[1] = R1 ;				\
  R0 = a[2] * b[0] + a[3] * b[2] ;			\
  R1 = a[2] * b[1] + a[3] * b[3] ;			\
  a[2] = R0 ; a[3] = R1 ;				
#else
  // slow and stupid loopy version
  int i , j , m ;
  GLU_complex R[ NC ] ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j ++ ) {
      R[j] = 0.0 ;
      for( m = 0 ; m < NC ; m ++  ) {
        R[ j ] += a[ m + i*NC ] * b[ j + m*NC ] ;
      }	
    }
    // copy back over to a ...
    for( j = 0 ; j < NC ; j ++ ) {
      a[ j + i*NC ] = R[j] ;
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
  // slow and stupid version
  size_t i , j , m ;
  register __m128d sum ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = _mm_setzero_pd( ) ;
      for( m = 0 ; m < NC ; m++  ) {
	sum = _mm_add_pd( sum , SSE2_MUL( *( B + m + NC * i ) , 
					  *( C + j + m * NC ) ) ) ;
      }
      A[ j + NC*i ] = sum ;
    }
  }
#endif
  return ;
}
#endif

#ifndef multab_suNC

void //__attribute__((hot))
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
  *( A + 2 ) = SSE2_FLIP( SSE2_CONJ( *( A + 1 ) ) ) ; 
  // a[3] =  conj( a[0] )
  *( A + 3 ) = SSE2_CONJ( *( A + 0 ) ) ;
#else
  return multab_suNC( a , b , c ) ;
#endif
  return ;
}

#endif

#endif // end of HAVE_IMMINTRIN_H

