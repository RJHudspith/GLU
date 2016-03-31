/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (invert.c) is part of GLU.

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
   @file invert.c
   @brief simple, arithmetic version of the inverse of a matrix using gauss-jordan
 */

#include "Mainfile.h"

// if we have Lapacke then we just let it do the heavy lifting
#ifdef HAVE_LAPACKE_H
  #include <lapacke.h> // we must be careful here
#endif 

//#define CLASSICAL_ADJOINT_INV  -> v. slow cofactor matrix version
//#define FULL_PIVOT             -> column and row pivoting

#if !(defined CLASSICAL_ADJOINT_INV) && ( NC > 2 ) 

// column elimination
static void
eliminate_column( double complex *__restrict a , 
		  double complex *__restrict inverse , 
		  const size_t i , 
		  const size_t j )
{
  const double complex fac1 = a[ i + j*NC ] / a[ i*(NC+1) ] ;  
  double complex *pA = a + i*NC ;
  size_t k ;
  // such a pattern elimintates cancelling zeros
  for( k = i + 1 ; k < NC ; k++ ) {
    a[ k + j*NC ] -= creal( fac1 ) * creal( pA[k] ) - cimag( fac1 ) * cimag( pA[k] ) 
      + I * ( creal( fac1 ) * cimag( pA[k] ) + cimag( fac1 ) * creal( pA[k] ) ) ;
  }
  pA = inverse + i*NC ;
  for( k = 0 ; k < NC ; k++ ) {
    // whatever we do to a, we do to the identity
    inverse[ k + j*NC ] -= creal( fac1 ) * creal( pA[k] ) - cimag( fac1 ) * cimag( pA[k] ) 
    + I * ( creal( fac1 ) * cimag( pA[k] ) + cimag( fac1 ) * creal( pA[k] ) ) ;
  }
  return ;
}

#ifdef FULL_PIVOT
// col swapper
static void
swap_cols( double complex *__restrict a , 
	   const size_t col_idx , 
	   const size_t piv )
{
  size_t l ;
  for( l = 0 ; l < NC ; l++ ) {
    register double complex tmp = a[col_idx+l*NC] ;
    a[col_idx+l*NC] = a[piv+l*NC] ;
    a[piv+l*NC] = tmp ;
  }
  return ;
}
#endif

// row swapper
static void
swap_rows( double complex *__restrict a , 
	   double complex *__restrict inverse , 
	   const size_t row_idx , 
	   const size_t piv )
{
  size_t l ;
  for( l = 0 ; l < NC ; l++ ) {
    register double complex tmp = a[l+row_idx*NC] ;
    a[l+row_idx*NC] = a[l+piv*NC] ;
    a[l+piv*NC] = tmp ;
    
    tmp = inverse[l+row_idx*NC] ;
    inverse[l+row_idx*NC] = inverse[l+piv*NC] ;
    inverse[l+piv*NC] = tmp ;
  }
  return ;
}

static int
gauss_jordan( GLU_complex M_1[ NCNC ] , 
	      const GLU_complex M[ NCNC ] )
{
  double complex a[ NCNC ] GLUalign , ainv[ NCNC ] GLUalign ;
  double row_best , col_best ;
  size_t i, j , col_piv ;

  // equate the necessary parts into double complex precision
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = (double complex)M[ i ] ;
    ainv[ i ] = 0.0 ;
  }
  for( i = 0 ; i < NC ; i++ ) {
    ainv[ i*(NC+1) ] = 1.0 ;
  }

#ifdef FULL_PIVOT
  size_t col_idx[ NC ] = {} , row_piv ;
#endif

  for( i = 0 ; i < NC-1 ; i++ ) {

    const size_t piv_idx = i*(NC+1) ;

    const double pivot = creal( a[piv_idx] ) * creal ( a[piv_idx] ) 
      +  cimag( a[piv_idx] ) * cimag( a[piv_idx] ) ;
    row_best = col_best = pivot ;
    col_piv = piv_idx ;

    #ifdef FULL_PIVOT
    row_piv = piv_idx ;
    col_idx[i] = i ;
    #endif
 
    size_t j ;
    for( j = i+1 ; j < NC ; j++ ) {
      #ifdef FULL_PIVOT
      const double row_attempt = creal( a[j+i*NC] ) * creal ( a[j+i*NC] ) 
	+  cimag( a[j+i*NC] ) * cimag( a[j+i*NC] ) ;
      //cabs( a[j+i*NC] ) ;
      if( row_attempt > row_best ) {
	row_best =  row_attempt ;
	row_piv = j ;
      }
      #endif
      const double col_attempt = creal( a[i+j*NC] ) * creal ( a[i+j*NC] ) 
	+  cimag( a[i+j*NC] ) * cimag( a[i+j*NC] ) ;
      if( col_attempt > col_best ) {
	col_best = col_attempt ;
	col_piv = j ;
      }
    }

    if( col_best > row_best && col_best > pivot ) {
      swap_rows( a , ainv , col_piv , i ) ;
    }
    #ifdef FULL_PIVOT
    if( col_best < row_best && row_best > pivot ) {
      //need to build up a permutation index here
      col_idx[i] = row_piv ;
      swap_cols( a , row_piv , i ) ;
    }
    #endif

    // perform gaussian elimination to obtain the upper triangular
    for( j = NC-1 ; j > i ; j-- ) { // go up in other columns
      eliminate_column( a , ainv , i , j ) ;
    }
  }

  // a is upper triangular, do the same for the upper half
  // no pivoting to be done here
  for( i = NC-1 ; i > 0 ; i-- ) {
    for( j = 0 ; j < i ; j++ ) {
      eliminate_column( a , ainv , i , j ) ;
    }
  }

  // multiply each row by its M_1 diagonal
  double complex am1 ;
  for( j = 0 ; j < NC ; j++ ) {
    if( creal( a[ j*( NC+1 ) ] ) == 0.0 && 
	cimag( a[ j*( NC+1 ) ] ) == 0.0 ) { 
      fprintf( stderr , "[INVERSE] Matrix is singular!!\n" ) ; 
      write_matrix( M ) ;
      return GLU_FAILURE ; 
    }
    am1 = 1.0 / a[ j*( NC+1 ) ] ;
    for( i = 0 ; i < NC ; i++ ) {
      ainv[ i + j*NC ] = creal( ainv[ i + j*NC ] ) * creal( am1 ) -
	cimag( ainv[ i + j*NC ] ) * cimag( am1 ) 
	+ I * ( creal( ainv[ i + j*NC ] ) * cimag( am1 ) + 
		cimag( ainv[ i + j*NC ] ) * creal( am1 )  ) ;
    }
  }

  #ifdef FULL_PIVOT
  // if we swap columns we must swap the rows of the answer
  for( i = NC-2 ; i >= 0 ; i-- ) {
    if( col_idx[i] != i ) {
      swap_rows( a , ainv , col_idx[i] , i ) ;
    }
  }
  #endif

  // push the double complex precision matrix back into whatever precision
  for( i = 0 ; i < NCNC ; i++ ) {
    M_1[i] = (GLU_complex)ainv[i] ;
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
  if( unlikely( cabs( deter ) < NC * PREC_TOL ) ) {
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
  M_1[ 0 ] = M[ 3 ] * INV_detM ;
  M_1[ 1 ] = -M[ 1 ] * INV_detM ;
  M_1[ 2 ] = -M[ 2 ] * INV_detM ;
  M_1[ 3 ] = M[ 0 ] * INV_detM ;
  return GLU_SUCCESS ;
  #else 
  return gauss_jordan( M_1 , M ) ;
  #endif
#endif 
}

// can use the Newton iteration as a cheaper alternative to actually computing the inverse
void
newton_approx_inverse( GLU_complex Zinv[ NCNC ] , // is Z_{k-1}^{-1}
		       const GLU_complex Z[ NCNC ] )
{
  GLU_complex BX[ NCNC ] GLUalign ;
  size_t iters , j ;
  for( iters = 0 ; iters < 4 ; iters++ ) {
    multab( BX , Z , Zinv ) ;
    // should be the identity
    if( fabs( creal( trace( BX ) ) - NC ) < PREC_TOL ) {
      break ;
    }
    multab_atomic_left( BX , Zinv ) ;
    for( j = 0 ; j < NCNC ; j++ ) {
      Zinv[ j ] = 2.0 * Zinv[ j ] - BX[ j ] ;
    }
  }
  return ;
}

// clear this one
#ifdef CLASSICAL_ADJOINT_INV
  #undef CLASSICAL_ADJOINT_INV
#endif

// clear the full pivoting, minimal improvement in using this
#ifdef FULL_PIVOT  
  #uncdf FULL_PIVOT
#endif
