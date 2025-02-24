/*
Copyright 2013-2025 Renwick James Hudspith

    This file (invert_SSE.c) is part of GLU.

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
   @file invert_SSE.c
   @brief gauss-jordan elimination with partial row pivoting, vectorised
 */
#include "Mainfile.h"

#if (defined HAVE_IMMINTRIN_H) && !(defined SINGLE_PREC)

// if we have Lapacke then we just let it do the heavy lifting
#ifdef HAVE_LAPACKE_H
  #include <lapacke.h> // we must be careful here
#endif 

#if !(defined CLASSICAL_ADJOINT_INV) && ( NC > 2 ) 

#include <immintrin.h>
#include "SSE2_OPS.h"

inline static __m128d
absfac( const __m128d a )
{
  register const __m128d z2 = _mm_shuffle_pd( a , a , 1 ) ;
  return _mm_add_pd( _mm_mul_pd( a , a ) , _mm_mul_pd( z2 , z2 ) ) ;
}

inline static __m128d
SSE2_inverse( const __m128d a )
{
  return _mm_div_pd( SSE2_CONJ( a ) , absfac( a ) ) ;
}

// column elimination
static void
eliminate_column( __m128d *a , 
		  __m128d *inverse ,
		  const __m128d fac ,
		  const size_t i , 
		  const size_t j )
{
  // common factor crops up everywhere
  register const __m128d fac1 = SSE2_MUL( a[ i + j*NC ] , fac ) ;
  __m128d *pA = ( a + i*(NC+1) + 1 ) , *ppA = ( a  + i + 1 + j*NC ) ;
  size_t k ;
  // such a pattern elimintates cancelling zeros
  for( k = i + 1 ; k < NC ; k++ ) {
    *ppA = _mm_sub_pd( *ppA , SSE2_MUL( fac1 , *pA ) ) ;
    pA ++ , ppA++ ;
  }
  // whatever we do to a, we do to the identity
  pA = ( inverse + i*NC ) ; ppA = ( inverse + j*NC ) ;
  for( k = 0 ; k < NC ; k++ ) {
    *ppA = _mm_sub_pd( *ppA , SSE2_MUL( fac1 , *pA ) ) ;
    pA ++ , ppA++ ;
  }
  return ;
}

// row swapper
static void
swap_rows( __m128d *a , 
	   __m128d *inverse , 
	   const size_t row_idx , 
	   const size_t piv )
{
  __m128d *pA = ( a + row_idx*NC ) , *ppA = ( a + piv*NC ) ;
  __m128d *pB = ( inverse + row_idx*NC ) , *ppB = ( inverse + piv*NC ) ;
  register __m128d tmp ;
  size_t l ;
  for( l = 0 ; l < NC ; l++ ) {
    tmp = *pA ; *pA = *ppA ; *ppA = tmp ;
    pA++ , ppA++ ;

    tmp = *pB ; *pB = *ppB ; *ppB = tmp ;
    pB++ , ppB++ ;
  }
  return ;
}

static int
gauss_jordan( GLU_complex M_1[ NCNC ] , 
	      const GLU_complex M[ NCNC ] )
{
  __m128d a[ NCNC ] GLUalign ; // temporary space to overwrite matrix
  register __m128d best , attempt , m1 , fac ;
  size_t i , j , piv ;

  // equate the necessary parts into double complex precision
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = _mm_setr_pd( creal( M[i] ) , cimag( M[i] ) ) ;
    M_1[ i ] = ( i%(NC+1) ) ? 0.0 :1.0 ;
  }

  // set these pointers, pB will be the inverse
  __m128d *pB = (__m128d*)M_1 , *pA = (__m128d*)a ;
  
  // loop over diagonal of the square matrix M
  for( i = 0 ; i < NC-1 ; i++ ) {

    // column pivot by selecting the largest in magnitude value
    piv = i ;
    best = absfac( *( pA + i*(NC+1) ) ) ;
    for( j = i+1 ; j < NC ; j++ ) {
       attempt = absfac( *( pB + i + j*NC ) ) ;
      if( _mm_ucomilt_sd( best , attempt ) ) { 
	piv = j ; 
	best = attempt ; 
      }
    }

    // if we must pivot then we do
    if( piv != i ) {
      swap_rows( pA , pB , piv , i ) ;
    }
  
    // perform gaussian elimination to obtain the upper triangular
    fac = _mm_div_pd( SSE2_CONJ( *( pA + i*(NC+1) ) ) , best ) ;
    for( j = NC-1 ; j > i ; j-- ) { // go up in other columns
      eliminate_column( pA , pB , fac , i , j ) ;
    }
  }

  // a is upper triangular, do the same for the upper half
  // no pivoting to be done here
  for( i = NC-1 ; i > 0 ; i-- ) {
    fac = SSE2_inverse( *( pA + i*(NC+1) ) ) ;
    for( j = 0 ; j < i ; j++ ) {
      eliminate_column( pA , pB , fac , i , j ) ;
    }
  }

  // multiply each row by its M_1 diagonal
  for( j = 0 ; j < NC ; j++ ) {
    m1 = SSE2_inverse( *pA ) ;
    for( i = 0 ; i < NC ; i++ ) {
      *pB = SSE2_MUL( *pB , m1 ) ;
      pB++ ;
    }
    pA += NC+1 ;
  }

  return GLU_SUCCESS ;
}
#endif

// calculates the inverse of the matrix M outputs to M_1 //
int 
inverse( GLU_complex M_1[ NCNC ] , 
	 const GLU_complex M[ NCNC ] )
{
#if (defined HAVE_LAPACKE_H)
  const int n = NC , lda = NC ;
  int ipiv[ NC + 1 ] ;
  memcpy( M_1 , M , NCNC * sizeof( GLU_complex ) ) ;
  int info = LAPACKE_zgetrf( LAPACK_ROW_MAJOR , n , n , M_1 , lda , ipiv ) ;
  info = LAPACKE_zgetri( LAPACK_ROW_MAJOR , n , M_1 , lda, ipiv ) ;
  return info ;
#elif (defined CLASSICAL_ADJOINT_INV)
  // define the adjunct //
  GLU_complex adjunct[ NCNC ] GLUalign ; 
  register GLU_complex deter = cofactor_transpose( adjunct , M ) ;
  // here we worry about numerical stability //
  if( cabs( deter ) < NC * PREC_TOL ) {
    fprintf( stderr , "[INVERSE] Matrix is singular !!! "
	     "deter=%1.14e %1.14e \n" , creal( deter ) , cimag( deter ) ) ; 
    write_matrix( M ) ; 
    return GLU_FAILURE ; 
  }
  // obtain inverse of M from 1/( detM )*adj( M ) //
  size_t i ;
  deter = 1.0 / deter ;
  for( i = 0 ; i < NCNC ; i++ ) { M_1[i] = adjunct[i] * deter ; }
  return GLU_SUCCESS ;
#else
  #if NC == 2 
  // use the identity, should warn for singular matrices
  const double complex INV_detM = 1.0 / ( det( M ) ) ;
  M_1[ 0 ] =  M[ 3 ] * INV_detM ;
  M_1[ 1 ] = -M[ 1 ] * INV_detM ;
  M_1[ 2 ] = -M[ 2 ] * INV_detM ;
  M_1[ 3 ] =  M[ 0 ] * INV_detM ;
  return GLU_SUCCESS ;
  #else 
  return gauss_jordan( M_1 , M ) ;
  #endif
#endif 
}

#endif // HAVE_IMMINTRIN_H
