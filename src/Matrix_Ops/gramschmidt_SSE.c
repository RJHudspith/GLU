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

   TODO :: speed up NC-generic gram-schmidt
 */

#include "Mainfile.h"

#include "GLU_rng.h"    // generate_NCxNC() is called

#if (defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC)

#include <immintrin.h>
#include "SSE2_OPS.h"

#if NC > 3
  #include "LU.h"         // LU_det()
#endif

//reunitarise an SU(3) matrix sped up for our requirements
void 
gram_reunit( GLU_complex *__restrict U )
{
#if NC == 3
    __m128d *u = ( __m128d* )U ;
  // orthogonalise the first row
  register __m128d sum ;

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
    __m128d *u = ( __m128d* )U ;

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
  size_t i , j , k , idx ;
  GLU_complex array[ (NC-1)*(NC-1) ] ;
  // perform the modified gram-schmidt
  register __m128d *v1 , *v2 ;
  register __m128d proj , norm , z1 ;
  for( k = 0 ; k < NC-1 ; k++ ) { // loop rows
    // do the projection
    for( i = 0 ; i < k ; i++ ) { // loop cols
      proj = _mm_setzero_pd( ) ;
      v1 = (__m128d*)U + k*NC ; v2 = (__m128d*)U + i*NC ;
      for( j = 0 ; j < NC ; j ++ ) {
	proj = _mm_add_pd( proj , SSE2_MUL_CONJ( *v1 , *v2 ) ) ;
	v1++ , v2++ ;
      }
      v1 = (__m128d*)U + k*NC ; v2 = (__m128d*)U + i*NC ;
      for( j = 0 ; j < NC ; j++ ) {
	*v1 = _mm_sub_pd( *v1 , SSE2_MUL( proj , *v2 ) ) ;
	v1++ , v2++ ;
      }
    }
    // norm code
    v1 = (__m128d*)U + k*NC ;
    norm = _mm_setzero_pd( ) ;
    for( j = 0 ; j < NC ; j++ ) {
      z1 = _mm_shuffle_pd( *v1 , *v1 , 1 ) ;
      norm = _mm_add_pd( norm , 
			 _mm_add_pd( _mm_mul_pd( *v1 , *v1 ) ,
				     _mm_mul_pd(  z1 ,  z1 ) ) ) ;
      v1++ ;
    }
    v1 = (__m128d*)U + k*NC ;
    norm = _mm_div_pd( _mm_setr_pd( 1.0 , 1.0 ) , _mm_sqrt_pd( norm ) ) ;
    for( j = 0 ; j < NC ; j++ ) {
      *v1 = _mm_mul_pd( *v1 , norm ) ; v1++ ;
    }
  }
  // initialise submatrix and compute sub-determinant
  for( k = 0 ; k < NC ; k++ ) {
    idx = 0 ;
    for( i = 0 ; i < NC-1 ; i++ ) {
      for( j = 0 ; j < NC ; j++ ) {
	if( j==k ) continue ;
	array[ idx++ ] = U[ j + i*NC ] ;
      }
    }
#if ( NC%2 == 0 )
    register const GLU_real mulfact = !(k&1) ? -1 :  1 ; 
#else
    register const GLU_real mulfact = !(k&1) ?  1 : -1 ; 
#endif
    // compute the determinant of the minors using the LU decomp
    U[ k + i*NC ] = conj( mulfact * (GLU_complex)LU_det_overwrite( NC-1 , array ) ) ;
  }
#endif
  return ;
}

// generate a random SU(N) matrix
void 
Sunitary_gen( GLU_complex Z[ NCNC ] )
{
  generate_NCxNC( Z ) ; // generate gaussian distributed elements of matrix
  gram_reunit( Z ) ;
  size_t counter = 0 ;
  while( is_unitary( Z ) == GLU_FALSE ) {
    fprintf( stderr , "[GRAM] not unitary! Redoing " ) ;
    generate_NCxNC( Z ) ; 
    gram_reunit( Z ) ;
    counter ++ ;
    // catastrophic failure
    if( counter > 10 ) exit(1) ;
  }
  return ;
}

#endif // <immintrin.h>
