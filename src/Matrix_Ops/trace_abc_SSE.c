/**
   @file trace_abc_SSE.c
   @brief trace of the product of three matrices (SSE2 versions)
 */
#include "Mainfile.h"

#if (defined HAVE_IMMINTRIN_H) && !(defined SINGLE_PREC)

#include <immintrin.h>
#include "SSE2_OPS.h"

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
  // last row is always a completion
  const __m128d a6 = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( a1 , a5 ) ,
					    SSE2_MUL( a2 , a4 ) ) ) ;
  const __m128d a7 = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( a2 , a3 ) ,
					    SSE2_MUL( a0 , a5 ) ) ) ;
  const __m128d a8 = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( a0 , a4 ) ,
					    SSE2_MUL( a1 , a3 ) ) ) ;
  // and compute the real part of the trace
  register __m128d sum ;
  sum = _mm_add_pd( _mm_mul_pd( a0 , *( C + 0 ) ) ,
		    _mm_add_pd( _mm_mul_pd( a1 , *( C + 1 ) ) ,
				_mm_mul_pd( a2 , *( C + 2 ) ) ) ) ;
  sum = _mm_add_pd( sum , _mm_add_pd( _mm_mul_pd( a3 , *( C + 3 ) ) ,
				      _mm_mul_pd( a4 , *( C + 4 ) ) ) ) ;
  sum = _mm_add_pd( sum , _mm_add_pd( _mm_mul_pd( a5 , *( C + 5 ) ) ,
				      _mm_mul_pd( a6 , *( C + 6 ) ) ) ) ;
  sum = _mm_add_pd( sum , _mm_add_pd( _mm_mul_pd( a7 , *( C + 7 ) ) ,
				      _mm_mul_pd( a8 , *( C + 8 ) ) ) ) ;
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
  //printf( "%1.15f\n" , creal( csum ) + cimag( csum ) ) ;
  //exit(1) ;
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
