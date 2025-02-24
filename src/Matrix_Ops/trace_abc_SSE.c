/*
Copyright 2013-2025 Renwick James Hudspith

    This file (trace_abc_SSE.c) is part of GLU.

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
   @file trace_abc_SSE.c
   @brief trace of the product of three matrices (SSE2 versions)
 */
#include "Mainfile.h"

#if (defined HAVE_IMMINTRIN_H) && !(defined SINGLE_PREC)

#include <immintrin.h>
#include "SSE2_OPS.h"

#define TEST

// compute the real part of the product of 3 su(3) matrices
double
Re_trace_abc_dag_suNC( const GLU_complex a[ NCNC ] , 
		       const GLU_complex b[ NCNC ] , 
		       const GLU_complex c[ NCNC ] )
{
  const __m128d *A = ( const __m128d* )a ;
  const __m128d *B = ( const __m128d* )b ;
  const __m128d *C = ( const __m128d* )c ;
  double complex csum ;
#if NC == 3
  register __m128d E , F , G , H , J , K ;
  E = SSE2_MUL( *( A + 0 ) , *( B + 0 ) ) ;
  F = SSE2_MUL( *( A + 1 ) , *( B + 3 ) ) ;
  G = SSE2_MUL( *( A + 2 ) , *( B + 6 ) ) ;
  H = SSE2_MUL( *( A + 0 ) , *( B + 1 ) ) ;
  J = SSE2_MUL( *( A + 1 ) , *( B + 4 ) ) ;
  K = SSE2_MUL( *( A + 2 ) , *( B + 7 ) ) ;
  E = _mm_add_pd( E , F ) ;
  H = _mm_add_pd( H , J ) ;
  const __m128d a0 = _mm_add_pd( E , G ) ;
  const __m128d a1 = _mm_add_pd( H , K ) ;

  E = SSE2_MUL( *( A + 0 ) , *( B + 2 ) ) ;
  F = SSE2_MUL( *( A + 1 ) , *( B + 5 ) ) ;
  G = SSE2_MUL( *( A + 2 ) , *( B + 8 ) ) ;
  H = SSE2_MUL( *( A + 3 ) , *( B + 0 ) ) ;
  J = SSE2_MUL( *( A + 4 ) , *( B + 3 ) ) ;
  K = SSE2_MUL( *( A + 5 ) , *( B + 6 ) ) ;
  E = _mm_add_pd( E , F ) ;
  H = _mm_add_pd( H , J ) ;
  const __m128d a2 = _mm_add_pd( E , G ) ;
  const __m128d a3 = _mm_add_pd( H , K ) ;

  E = SSE2_MUL( *( A + 3 ) , *( B + 1 ) ) ;
  F = SSE2_MUL( *( A + 4 ) , *( B + 4 ) ) ;
  G = SSE2_MUL( *( A + 5 ) , *( B + 7 ) ) ;
  H = SSE2_MUL( *( A + 3 ) , *( B + 2 ) ) ;
  J = SSE2_MUL( *( A + 4 ) , *( B + 5 ) ) ;
  K = SSE2_MUL( *( A + 5 ) , *( B + 8 ) ) ;
  E = _mm_add_pd( E , F ) ;
  H = _mm_add_pd( H , J ) ;
  const __m128d a4 = _mm_add_pd( E , G ) ;
  const __m128d a5 = _mm_add_pd( H , K ) ;

  // this does the completion
  F = SSE2_MUL( a1 , a3 ) ;
  H = SSE2_MUL( a0 , a4 ) ;
  G = SSE2_MUL( a2 , a4 ) ;
  K = SSE2_MUL( a1 , a5 ) ;
  E = SSE2_MUL( a2 , a3 ) ;
  J = SSE2_MUL( a0 , a5 ) ;
  H = _mm_sub_pd( H , F ) ;
  K = _mm_sub_pd( K , G ) ;
  E = _mm_sub_pd( E , J ) ;
  const __m128d a8 = SSE2_CONJ( H ) ;
  const __m128d a6 = SSE2_CONJ( K ) ;
  const __m128d a7 = SSE2_CONJ( E ) ;
    
  // and compute the real part of the trace
  E =  _mm_mul_pd( a2 , *( C + 2 ) ) ;
  F =  _mm_mul_pd( a1 , *( C + 1 ) ) ;
  G =  _mm_mul_pd( a0 , *( C + 0 ) ) ;
  E = _mm_add_pd( E , F ) ;
  E = _mm_add_pd( E , G ) ;

  register __m128d sum = _mm_setzero_pd( ) ;
  sum = _mm_add_pd( sum , E ) ;

  H =  _mm_mul_pd( a3 , *( C + 3 ) ) ;
  J =  _mm_mul_pd( a4 , *( C + 4 ) ) ;
  K =  _mm_mul_pd( a5 , *( C + 5 ) ) ;
  H = _mm_add_pd( H , J ) ;
  H = _mm_add_pd( H , K ) ;  
  sum = _mm_add_pd( sum , H ) ;

  E =  _mm_mul_pd( a6 , *( C + 6 ) ) ;
  F =  _mm_mul_pd( a7 , *( C + 7 ) ) ;
  G =  _mm_mul_pd( a8 , *( C + 8 ) ) ;
  E = _mm_add_pd( E , F ) ;
  E = _mm_add_pd( E , G ) ;
  sum = _mm_add_pd( sum , E ) ;
  
  _mm_store_pd( (void*)&csum , sum ) ;
#elif NC == 2
  // puts the four parts of the sum into the upper and lower parts of 
  // two SSE-packed double words
  register const __m128d a0 = _mm_sub_pd( SSE2_MUL( *( A + 0 ) , *( B + 0 ) ) ,
					  SSE2_MUL_CONJ( *( A + 1 ) , *( B + 1 ) ) ) ;
  register const __m128d a1 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 1 ) ) ,
					  SSE2_MUL_CONJ( *( A + 1 ) , *( B + 0 ) ) ) ;
  register const __m128d sum = _mm_add_pd( _mm_mul_pd( a0 , *( C + 0 ) ) , 
					   _mm_mul_pd( a1 , *( C + 1 ) ) ) ;
  // multiply sum by 2
  _mm_store_pd( (void*)&csum , _mm_add_pd( sum , sum ) ) ; 
#else
  register __m128d sum = _mm_setzero_pd( ) , insum ;
  size_t i , j , k ;
  for( i = 0 ; i < NC ; i++ ) {
    A = (const __m128d*)a ;
    for( j = 0 ; j < NC ; j++ ) {
      B = (const __m128d*)b ;
      insum = _mm_setzero_pd( ) ;
      for( k = 0 ; k < NC ; k++ ) {
	insum = _mm_add_pd( insum , SSE2_MUL( A[k] , B[i] ) ) ;
	B += NC ;
      }
      sum = _mm_add_pd( sum , _mm_mul_pd( insum , C[i+j*NC] ) ) ;
      // increment our pointers
      A += NC ;
    }
  }
  // multiply by 2
  _mm_store_pd( (void*)&csum , sum ) ;
#endif
  return creal( csum ) + cimag( csum ) ;  
}

// Trace of the product of three matrices //
void
trace_abc( GLU_complex *__restrict tr , 
	   const GLU_complex a[ NCNC ] , 
	   const GLU_complex b[ NCNC ] , 
	   const GLU_complex c[ NCNC ] )
{
  const __m128d *A = ( const __m128d* )a ;
  const __m128d *B = ( const __m128d* )b ;
  const __m128d *C = ( const __m128d* )c ;
  double complex csum ;
#if NC == 3
  // first row
  const __m128d a0 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 0 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 1 ) , *( B + 3 ) ) ,
					     SSE2_MUL( *( A + 2 ) , *( B + 6 ) ) ) ) ;
  const __m128d a1 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 1 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 1 ) , *( B + 4 ) ) ,
					     SSE2_MUL( *( A + 2 ) , *( B + 7 ) ) ) ) ;
  const __m128d a2 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 2 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 1 ) , *( B + 5 ) ) ,
					     SSE2_MUL( *( A + 2 ) , *( B + 8 ) ) ) ) ;
  // second row
  const __m128d a3 = _mm_add_pd( SSE2_MUL( *( A + 3 ) , *( B + 0 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 4 ) , *( B + 3 ) ) ,
					     SSE2_MUL( *( A + 5 ) , *( B + 6 ) ) ) ) ;
  const __m128d a4 = _mm_add_pd( SSE2_MUL( *( A + 3 ) , *( B + 1 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 4 ) , *( B + 4 ) ) ,
					     SSE2_MUL( *( A + 5 ) , *( B + 7 ) ) ) ) ;
  const __m128d a5 = _mm_add_pd( SSE2_MUL( *( A + 3 ) , *( B + 2 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 4 ) , *( B + 5 ) ) ,
					     SSE2_MUL( *( A + 5 ) , *( B + 8 ) ) ) ) ;
  // third row
  const __m128d a6 = _mm_add_pd( SSE2_MUL( *( A + 6 ) , *( B + 0 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 7 ) , *( B + 3 ) ) ,
					     SSE2_MUL( *( A + 8 ) , *( B + 6 ) ) ) ) ;
  const __m128d a7 = _mm_add_pd( SSE2_MUL( *( A + 6 ) , *( B + 1 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 7 ) , *( B + 4 ) ) ,
					     SSE2_MUL( *( A + 8 ) , *( B + 7 ) ) ) ) ;
  const __m128d a8 = _mm_add_pd( SSE2_MUL( *( A + 6 ) , *( B + 2 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 7 ) , *( B + 5 ) ) ,
					     SSE2_MUL( *( A + 8 ) , *( B + 8 ) ) ) ) ;
  // and compute the real part of the trace
  register __m128d sum = _mm_setzero_pd( ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a0 , *( C + 0 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a1 , *( C + 3 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a2 , *( C + 6 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a3 , *( C + 1 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a4 , *( C + 4 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a5 , *( C + 7 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a6 , *( C + 2 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a7 , *( C + 5 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a8 , *( C + 8 ) ) ) ;
#elif NC == 2
  // first row of the multiply
  const __m128d a0 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 0 ) ) ,
				 SSE2_MUL( *( A + 1 ) , *( B + 2 ) ) ) ;
  const __m128d a1 = _mm_add_pd( SSE2_MUL( *( A + 2 ) , *( B + 0 ) ) ,
				 SSE2_MUL( *( A + 3 ) , *( B + 2 ) ) ) ;
  // second row of the multiply
  const __m128d a2 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 1 ) ) ,
				 SSE2_MUL( *( A + 1 ) , *( B + 3 ) ) ) ;
  const __m128d a3 = _mm_add_pd( SSE2_MUL( *( A + 2 ) , *( B + 1 ) ) ,
				 SSE2_MUL( *( A + 3 ) , *( B + 3 ) ) ) ;
  // compute the sum
  register __m128d sum = _mm_setzero_pd( ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a0 , *( C + 0 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a1 , *( C + 1 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a2 , *( C + 2 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL( a3 , *( C + 3 ) ) ) ;
#else
  register __m128d sum = _mm_setzero_pd( ) ;
  size_t i , j , k ;
  for( i = 0 ; i < NC ; i++ ) {
    A = (const __m128d*)a ;
    for( j = 0 ; j < NC ; j++ ) {
      B = (const __m128d*)b ;
      register __m128d insum = _mm_setzero_pd( ) ;
      for( k = 0 ; k < NC ; k++ ) {
	insum = _mm_add_pd( insum , SSE2_MUL( A[k] , B[i] ) ) ;
	B += NC ;
      }
      sum = _mm_add_pd( sum , SSE2_MUL( insum , *C ) ) ;
      // increment our pointers
      A += NC , C++ ;
    }
  }
#endif
  _mm_store_pd( (void*)&csum , sum ) ;
  *tr = csum ;
  return ;
}

// this is trace( a . b . c^{\dagger} )
void
trace_abc_dag( GLU_complex *__restrict tr , 
	       const GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC ] , 
	       const GLU_complex c[ NCNC ] )
{
  const __m128d *A = ( const __m128d* )a ;
  const __m128d *B = ( const __m128d* )b ;
  const __m128d *C = ( const __m128d* )c ;
  double complex csum ;
#if NC == 3
  //const GLU_complex a0 = ( a[0] * b[0] + a[1] * b[3] + a[2] * b[6] ) ;
  const __m128d a0 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 0 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 1 ) , *( B + 3 ) ) ,
					     SSE2_MUL( *( A + 2 ) , *( B + 6 ) ) ) ) ;
  const __m128d a1 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 1 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 1 ) , *( B + 4 ) ) ,
					     SSE2_MUL( *( A + 2 ) , *( B + 7 ) ) ) ) ;
  const __m128d a2 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 2 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 1 ) , *( B + 5 ) ) ,
					     SSE2_MUL( *( A + 2 ) , *( B + 8 ) ) ) ) ;
  // second row
  const __m128d a3 = _mm_add_pd( SSE2_MUL( *( A + 3 ) , *( B + 0 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 4 ) , *( B + 3 ) ) ,
					     SSE2_MUL( *( A + 5 ) , *( B + 6 ) ) ) ) ;
  const __m128d a4 = _mm_add_pd( SSE2_MUL( *( A + 3 ) , *( B + 1 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 4 ) , *( B + 4 ) ) ,
					     SSE2_MUL( *( A + 5 ) , *( B + 7 ) ) ) ) ;
  const __m128d a5 = _mm_add_pd( SSE2_MUL( *( A + 3 ) , *( B + 2 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 4 ) , *( B + 5 ) ) ,
					     SSE2_MUL( *( A + 5 ) , *( B + 8 ) ) ) ) ;
  // third row
  const __m128d a6 = _mm_add_pd( SSE2_MUL( *( A + 6 ) , *( B + 0 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 7 ) , *( B + 3 ) ) ,
					     SSE2_MUL( *( A + 8 ) , *( B + 6 ) ) ) ) ;
  const __m128d a7 = _mm_add_pd( SSE2_MUL( *( A + 6 ) , *( B + 1 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 7 ) , *( B + 4 ) ) ,
					     SSE2_MUL( *( A + 8 ) , *( B + 7 ) ) ) ) ;
  const __m128d a8 = _mm_add_pd( SSE2_MUL( *( A + 6 ) , *( B + 2 ) ) ,
				 _mm_add_pd( SSE2_MUL( *( A + 7 ) , *( B + 5 ) ) ,
					     SSE2_MUL( *( A + 8 ) , *( B + 8 ) ) ) ) ;
  // and compute the real part of the trace
  register __m128d sum = _mm_setzero_pd( ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a0 , *( C + 0 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a1 , *( C + 1 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a2 , *( C + 2 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a3 , *( C + 3 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a4 , *( C + 4 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a5 , *( C + 5 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a6 , *( C + 6 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a7 , *( C + 7 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a8 , *( C + 8 ) ) ) ;
#elif NC == 2
  // first row of the multiply
  const __m128d a0 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 0 ) ) ,
				 SSE2_MUL( *( A + 1 ) , *( B + 2 ) ) ) ;
  const __m128d a1 = _mm_add_pd( SSE2_MUL( *( A + 0 ) , *( B + 1 ) ) ,
				 SSE2_MUL( *( A + 1 ) , *( B + 3 ) ) ) ;
  // second row of the multiply
  const __m128d a2 = _mm_add_pd( SSE2_MUL( *( A + 2 ) , *( B + 0 ) ) ,
				 SSE2_MUL( *( A + 3 ) , *( B + 2 ) ) ) ;
  const __m128d a3 = _mm_add_pd( SSE2_MUL( *( A + 2 ) , *( B + 1 ) ) ,
				 SSE2_MUL( *( A + 3 ) , *( B + 3 ) ) ) ;
  // compute the sum
  register __m128d sum = _mm_setzero_pd( ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a0 , *( C + 0 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a1 , *( C + 1 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a2 , *( C + 2 ) ) ) ;
  sum = _mm_add_pd( sum , SSE2_MUL_CONJ( a3 , *( C + 3 ) ) ) ;
#else
  register __m128d sum = _mm_setzero_pd( ) ;
  size_t i , j , k ;
  for( i = 0 ; i < NC ; i++ ) {
    A = (const __m128d*)a ;
    for( j = 0 ; j < NC ; j++ ) {
      B = (const __m128d*)b ;
      register __m128d insum = _mm_setzero_pd( ) ;
      for( k = 0 ; k < NC ; k++ ) {
	insum = _mm_add_pd( insum , SSE2_MUL( A[k] , B[i] ) ) ;
	//insum = _mm_add_pd( insum , SSE2_MUL( A[k] , B[j] ) ) ;
	B += NC ;
      }
      //sum = _mm_add_pd( sum , SSE2_MUL_CONJ( insum , *C ) ) ;
      sum = _mm_add_pd( sum , SSE2_MUL_CONJ( insum , C[i+j*NC] ) ) ;
      // increment our pointers
      A += NC ; //, C++ ;
    }
  }
#endif
  _mm_store_pd( (void*)&csum , sum ) ;
  *tr = (GLU_complex)csum ;
  return ;
}

// this is trace( a . b . c^{\dagger} )
void
trace_abc_dag_Re( GLU_real *__restrict tr , 
		  const GLU_complex a[ NCNC ] , 
		  const GLU_complex b[ NCNC ] , 
		  const GLU_complex c[ NCNC ] )
{
  GLU_complex Ctr ;
  trace_abc_dag( &Ctr , a , b , c ) ;
  *tr = creal( Ctr ) ;
  return ;
}

#endif
