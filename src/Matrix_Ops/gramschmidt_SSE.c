/*
    Copyright 2013 Renwick James Hudspith

    This file (gramschmidt_SSE.c) is part of GLU.

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
   @file gramschmidt.c
   @brief reunitarisation procedures SSEd for SU(3) and SU(2)
 */

#include "Mainfile.h"

#include "GLU_rng.h"  // generate_NCxNC() is called

#if (defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC)

#include <immintrin.h>
#include "SSE2_OPS.h"

#if NC > 3

// gramschmidt projection V = V - V.U^{\dagger}
static void
project( V , U )
     GLU_complex *__restrict V ;
     const GLU_complex *__restrict U ;
{
  int i ;
  register double projRE = 0.0 , projIM = 0.0 ;
  for( i = 0 ; i < NC ; i ++ ) {
    projRE += ( double )( creal( V[i] ) * creal( U[i] ) ) ;
    projRE += ( double )( cimag( V[i] ) * cimag( U[i] ) ) ;
    projIM += ( double )( creal( U[i] ) * cimag( V[i] ) ) ;
    projIM -= ( double )( cimag( U[i] ) * creal( V[i] ) ) ; 
  } 
  for( i = 0 ; i < NC ; i++ ) {
    V[i] -= projRE * creal( U[i] ) - projIM * cimag( U[i] ) +
      I * ( projRE * cimag( U[i] ) + projIM * creal( U[i] ) ) ;
  }
  return ;
}

// normalize a vector //
INLINE_STATIC_VOID
vect_norm2( a )
     GLU_complex *__restrict a ;
{
  int mu ;
  register double norm = 0.0 ;
  for( mu = 0 ; mu < NC ; mu++ ) {
    norm += (double)( creal( a[mu] ) * creal( a[mu] ) ) +	\
            (double)( cimag( a[mu] ) * cimag( a[mu] ) ) ;
  }
  norm = 1. / sqrt( norm ) ;
  for( mu = 0 ; mu < NC ; mu++ ) {
    a[ mu ] *= norm ;
  }
  return ;
}
#endif

//reunitarise an SU(3) matrix sped up for our requirements
void 
reunit2( GLU_complex *__restrict U)
{
    __m128d *u = ( __m128d* )U ;
#if NC == 3
  // orthogonalise the first row
  register __m128d sum ;

  double complex s ;

  // orthogonalise first row

  sum = _mm_add_pd( _mm_mul_pd( *( u + 0 ) , *( u + 0 ) ) ,
		    _mm_add_pd( _mm_mul_pd( *( u + 3 ) , *( u + 3 ) ) ,
				_mm_mul_pd( *( u + 6 ) , *( u + 6 ) ) ) ) ;

  sum = _mm_add_pd( sum , _mm_shuffle_pd( sum , sum , 1 ) ) ;
  sum = _mm_div_pd( _mm_setr_pd( 1. , 1. ) , _mm_sqrt_pd( sum ) ) ; // has the norm
 
  *( u + 0 ) = _mm_mul_pd( *( u + 0 ) , sum ) ;
  *( u + 3 ) = _mm_mul_pd( *( u + 3 ) , sum ) ;
  *( u + 6 ) = _mm_mul_pd( *( u + 6 ) , sum ) ;

  // gramschmidt second with respect to first
  // u = u - vw / ww
  sum = _mm_add_pd( SSE2_MUL_CONJ( *( u + 1 ) , *( u + 0 ) ) ,
		    _mm_add_pd( SSE2_MUL_CONJ( *( u + 4 ) , *( u + 3 ) ) ,
				SSE2_MUL_CONJ( *( u + 7 ) , *( u + 6 ) ) ) ) ;
  *( u + 1 ) = _mm_sub_pd( *( u + 1 ) , SSE2_MUL( sum , *( u + 0 ) ) ) ;
  *( u + 4 ) = _mm_sub_pd( *( u + 4 ) , SSE2_MUL( sum , *( u + 3 ) ) ) ;
  *( u + 7 ) = _mm_sub_pd( *( u + 7 ) , SSE2_MUL( sum , *( u + 6 ) ) ) ;

  // normalise
  sum = _mm_add_pd( _mm_mul_pd( *( u + 1 ) , *( u + 1 ) ) ,
		    _mm_add_pd( _mm_mul_pd( *( u + 4 ) , *( u + 4 ) ) ,
				_mm_mul_pd( *( u + 7 ) , *( u + 7 ) ) ) ) ;

  sum = _mm_add_pd( sum , _mm_shuffle_pd( sum , sum , 1 ) ) ;
  sum = _mm_div_pd( _mm_setr_pd( 1. , 1. ) , _mm_sqrt_pd( sum ) ) ; // has the norm

  *( u + 1 ) = _mm_mul_pd( *( u + 1 ) , sum ) ;
  *( u + 4 ) = _mm_mul_pd( *( u + 4 ) , sum ) ;
  *( u + 7 ) = _mm_mul_pd( *( u + 7 ) , sum ) ;

  // and complete using the minors
  *( u + 2 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( u + 3 ) , *( u + 7 ) ) ,
				      SSE2_MUL( *( u + 4 ) , *( u + 6 ) ) ) ) ;
  *( u + 5 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( u + 1 ) , *( u + 6 ) ) ,
				      SSE2_MUL( *( u + 0 ) , *( u + 7 ) ) ) ) ;
  *( u + 8 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( u + 0 ) , *( u + 4 ) ) ,
				      SSE2_MUL( *( u + 1 ) , *( u + 3 ) ) ) ) ;
#elif NC == 2

  /*
  GLU_real *uu = (GLU_real*)U ;
  const double v00R = *( uu + 0 ) ;
  const double v00I = *( uu + 1 ) ;
  const double v01R = *( uu + 4 ) ;
  const double v01I = *( uu + 5 ) ;

  double norm = sqrt( v00R * v00R + v00I * v00I + v01R * v01R + v01I * v01I ) ;
  norm = 1. / norm ;

  U[0] =  ( v00R * norm + I * v00I * norm ) ; 
  U[1] = -( v01R * norm - I * v01I * norm ) ; 
  U[2] =  ( v01R * norm + I * v01I * norm ) ; 
  U[3] =  ( v00R * norm - I * v00I * norm ) ; 
  */

  // compute the norm 
  register __m128d sum ;
  sum = _mm_add_pd( _mm_mul_pd( *( u + 0 ) , *( u + 0 ) ) ,
		    _mm_mul_pd( *( u + 2 ) , *( u + 2 ) ) ) ;
  sum = _mm_add_pd( sum , _mm_shuffle_pd( sum , sum , 1 ) ) ;
  sum = _mm_div_pd( _mm_setr_pd( 1 , 1 ) , _mm_sqrt_pd( sum ) ) ;
  *( u + 0 ) = _mm_mul_pd( *( u + 0 ) , sum ) ;
  *( u + 2 ) = _mm_mul_pd( *( u + 2 ) , sum ) ;
  *( u + 1 ) = SSE_FLIP( SSE2_CONJ( *( u + 2 ) ) ) ; 
  *( u + 3 ) = SSE2_CONJ( *( u + 0 ) ) ;

#else
  // need a modified gram-schmidt process, leaves the
  // bottom row untouched
  int i , j ;
  GLU_complex v[ NC - 1 ][ NC ] ;

  // equate our vector to our matrix , parallelism here is a bit extreme
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      v[i][j] = U[ j + ( NC * i ) ] ;
    }
  }

  // perform the modified gram schmidt
  int k ;
  for( k = 0 ; k < NC-1 ; k++ ) {
    for( i = 0 ; i < k ; i++ ) { 
      project( v[k] , v[i] ) ;      
    }
    // normalize
    vect_norm2( v[k] ) ;
    // and put back into U, need the index i to make sure the vector gets normalised
    for( j = 0 ; j < NC ; j++ ) {
      U[ j + ( NC * i ) ] = v[i][j] ;
    }
  }

  // OK so this process leaves the bottom row as was, need to complete
  GLU_complex *array = malloc( ( NC - 1 ) * ( NC - 1 ) * sizeof( GLU_complex ) ) ;
  for( i = (NCNC-NC) ; i < NCNC ; i++ ) { // our bona-fide minor index loops the bottom row
    int idx = 0 ;
    for( j = 0 ; j < ( NCNC - NC ) ; j++ ) {
      if( ( j%NC != i%NC ) ) { // remove columns  ( j/NC != i/NC ) is implicit!!
	// pack array
	array[idx++] = U[j] ;
      } 
    }
    // compute the determinant
    #if ( NC%2 == 0 )
    register const GLU_real mulfact = ( i % 2 == 0 ) ? -1.0 : 1.0 ; 
    #else
    register const GLU_real mulfact = ( i % 2 == 0 ) ? 1.0 : -1.0 ; 
    #endif
    // compute the determinant of the minors using the LU decomp
    U[i] = conj( mulfact * (GLU_complex)LU_det( NC-1 , array ) ) ;
  }
  free( array ) ;
#endif
  return ;
}

//reunitarise an SU(Nc) matrix
void 
reunit( GLU_complex Z[ NCNC ] ,
	const GLU_complex U[ NCNC ] ) 
{
  equiv( Z , U ) ;
  reunit2( Z ) ;
  return ;
}

// generate a random SU(N) matrix
void 
Sunitary_gen( GLU_complex Z[ NCNC ] )
{
  generate_NCxNC( Z ) ; // generate gaussian distributed elements of matrix
  reunit2( Z ) ;
  while( !is_unitary( Z ) ) {
    printf("not unitary! Redoing ") ;
    generate_NCxNC( Z ) ; 
    reunit2( Z ) ;
  }
  return ;
}

#endif // <immintrin.h>

