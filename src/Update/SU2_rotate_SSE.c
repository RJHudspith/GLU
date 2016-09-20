/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (SU2_rotate_SSE.h) is part of GLU.

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
   @file SU2_rotate_SSE.c
   @brief SU2 subgroup rotations on SU(#NC) matrices
 */
#include "Mainfile.h"

// have to have intrinsics turned on and be in double precision
#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC )

#include <immintrin.h>
#include "SSE2_OPS.h"

// compute the relevant su2 indices IIRC shamelessly stolen from chroma
void
compute_pertinent_indices( void )
{
  Latt.su2_data = malloc( NSU2SUBGROUPS * sizeof( struct su2_subgroups ) ) ;
  size_t su2_index , i1 , i2 ;
  for( su2_index = 0 ; su2_index < NSU2SUBGROUPS ; su2_index ++ ) {
    size_t found = 0 , del_i = 0 , index = 0 ;
    while ( del_i < (NC-1) && found == 0 ) {
      del_i++;
      for ( i1 = 0; i1 < (NC-del_i); i1++ ) {
	if ( index == su2_index ) {
	  found = 1;
	  break;
	}
	index++;
      }
    }
    i2 = i1 + del_i ;
    Latt.su2_data[ su2_index ].idx_a = i1 + NC * i1 ;
    Latt.su2_data[ su2_index ].idx_b = i2 + NC * i1 ;
    Latt.su2_data[ su2_index ].idx_c = i1 + NC * i2 ;
    Latt.su2_data[ su2_index ].idx_d = i2 + NC * i2 ;
  }
  return ;
}

// free the struct
void
free_su2_data( void )
{
  free( Latt.su2_data ) ;
}

// only compute the necessary indices of su2_i = subgroup( U*staple^\dagger )
void
only_subgroup( GLU_complex *s0 ,
	       GLU_complex *s1 ,
	       double *scale ,
	       const GLU_complex U[ NCNC ] ,
	       const GLU_complex staple[ NCNC ] ,
	       const size_t su2_index )
{
  const __m128d *u = (const __m128d*)U ;
  const __m128d *s = (const __m128d*)staple ;

  register __m128d sm0 ; 
  register __m128d sm1 ;
#if NC == 3
  switch( su2_index%3 ) { // I don't like this
    // rotation 1
    //  |  s0   s1  0 |
    //  | -s1*  s0* 0 |
    //  |  0     0  1 |
  case 0 :
    sm0 = _mm_add_pd(
		     // temp0
		     _mm_add_pd( SSE2_MUL_CONJ( *( u + 0 ) , *( s + 0 ) ) ,
				 _mm_add_pd( SSE2_MUL_CONJ( *( u + 1 ) , *( s + 1 ) ) ,
					     SSE2_MUL_CONJ( *( u + 2 ) , *( s + 2 ) ) ) ) ,
		      // temp3^*
		     _mm_add_pd( SSE2_MULCONJ( *( u + 3 ) , *( s + 3 ) ) ,
				 _mm_add_pd( SSE2_MULCONJ( *( u + 4 ) , *( s + 4 ) ) ,
					     SSE2_MULCONJ( *( u + 5 ) , *( s + 5 ) ) ) ) ) ;
    sm1 = _mm_sub_pd( 
		     // temp1
		     _mm_add_pd( SSE2_MUL_CONJ( *( u + 0 ) , *( s + 3 ) ) ,
				 _mm_add_pd( SSE2_MUL_CONJ( *( u + 1 ) , *( s + 4 ) ) ,
					     SSE2_MUL_CONJ( *( u + 2 ) , *( s + 5 ) ) ) ) ,
		     // temp2^*
		     _mm_add_pd( SSE2_MULCONJ( *( u + 3 ) , *( s + 0 ) ) ,
				 _mm_add_pd( SSE2_MULCONJ( *( u + 4 ) , *( s + 1 ) ) ,
					     SSE2_MULCONJ( *( u + 5 ) , *( s + 2 ) ) ) ) ) ;
    break ;
  case 1 :
    // rotation 2
    //  |  1    0   0  |
    //  |  0   s0  s1  |
    //  |  0  -s1* s0* |
    sm0 = _mm_add_pd( 
		     // temp0
		     _mm_add_pd( SSE2_MUL_CONJ( *( u + 3 ) , *( s + 3 ) ) ,
				 _mm_add_pd( SSE2_MUL_CONJ( *( u + 4 ) , *( s + 4 ) ) ,
					     SSE2_MUL_CONJ( *( u + 5 ) , *( s + 5 ) ) ) ) ,
		     // temp3^*
		     _mm_add_pd( SSE2_MULCONJ( *( u + 6 ) , *( s + 6 ) ) ,
				 _mm_add_pd( SSE2_MULCONJ( *( u + 7 ) , *( s + 7 ) ) ,
					     SSE2_MULCONJ( *( u + 8 ) , *( s + 8 ) ) ) ) ) ;
    sm1 = _mm_sub_pd(
		     // temp1
		     _mm_add_pd( SSE2_MUL_CONJ( *( u + 3 ) , *( s + 6 ) ) ,
				 _mm_add_pd( SSE2_MUL_CONJ( *( u + 4 ) , *( s + 7 ) ) ,
					     SSE2_MUL_CONJ( *( u + 5 ) , *( s + 8 ) ) ) ) ,
		     // temp2^*
		     _mm_add_pd( SSE2_MULCONJ( *( u + 6 ) , *( s + 3 ) ) ,
				 _mm_add_pd( SSE2_MULCONJ( *( u + 7 ) , *( s + 4 ) ) ,
					     SSE2_MULCONJ( *( u + 8 ) , *( s + 5 ) ) ) ) ) ;
    break ;
  case 2 :
    // rotation 3
    //  | s0*  0  -s1 |
    //  |  0   1   0  |
    //  | s1   0  s0  |
    sm0 = _mm_add_pd( 
		     // temp3^*
		     _mm_add_pd( SSE2_MULCONJ( *( u + 0 ) , *( s + 0 ) ) ,
				 _mm_add_pd( SSE2_MULCONJ( *( u + 1 ) , *( s + 1 ) ) ,
					     SSE2_MULCONJ( *( u + 2 ) , *( s + 2 ) ) ) ) ,
		     // temp0
		     _mm_add_pd( SSE2_MUL_CONJ( *( u + 6 ) , *( s + 6 ) ) ,
				 _mm_add_pd( SSE2_MUL_CONJ( *( u + 7 ) , *( s + 7 ) ) ,
					     SSE2_MUL_CONJ( *( u + 8 ) , *( s + 8 ) ) ) ) ) ;
    sm1 = _mm_sub_pd(
		     // temp1
		     _mm_add_pd( SSE2_MUL_CONJ( *( u + 6 ) , *( s + 0 ) ) ,
				 _mm_add_pd( SSE2_MUL_CONJ( *( u + 7 ) , *( s + 1 ) ) ,
					     SSE2_MUL_CONJ( *( u + 8 ) , *( s + 2 ) ) ) ) ,
		     // temp2^*
		     _mm_add_pd( SSE2_MULCONJ( *( u + 0 ) , *( s + 6 ) ) ,
				 _mm_add_pd( SSE2_MULCONJ( *( u + 1 ) , *( s + 7 ) ) ,
					     SSE2_MULCONJ( *( u + 2 ) , *( s + 8 ) ) ) ) ) ;
    break ;
  }
#elif NC == 2
  sm0 = _mm_add_pd( 
		   // temp0
		   _mm_add_pd( SSE2_MUL_CONJ( *( u + 0 ) , *( s + 0 ) ) , 
			       SSE2_MUL_CONJ( *( u + 1 ) , *( s + 1 ) ) ) ,
		   // temp3^*
		   _mm_add_pd( SSE2_MULCONJ( *( u + 2 ) , *( s + 2 ) ) , 
			       SSE2_MULCONJ( *( u + 3 ) , *( s + 3 ) ) ) ) ;
  
  sm1 = _mm_sub_pd( 
		   // temp1
		   _mm_add_pd( SSE2_MUL_CONJ( *( u + 0 ) , *( s + 2 ) ) , 
			       SSE2_MUL_CONJ( *( u + 1 ) , *( s + 3 ) ) ) ,
		   // temp2^*
		   _mm_add_pd( SSE2_MULCONJ( *( u + 2 ) , *( s + 0 ) ) ,
			       SSE2_MULCONJ( *( u + 3 ) , *( s + 1 ) ) ) ) ;
#else
  // su(N) version
  const size_t row_a = Latt.su2_data[ su2_index ].idx_a / NC ;
  const size_t col_b = Latt.su2_data[ su2_index ].idx_b % NC ;

  // prefetch the staple & link indices
  const __m128d *S1 = ( s + NC * row_a ) , *S2 = ( s + NC * col_b ) ; 
  const __m128d *U1 = ( u + NC * row_a ) , *U2 = ( u + NC * col_b ) ; 

  // initialise to zero & perform multiplication
  sm0 = _mm_setzero_pd() ; sm1 = _mm_setzero_pd() ;

  size_t i ;
  for( i = 0 ; i < NC ; i++ ) {
    sm0 = _mm_add_pd( sm0 ,
		      _mm_add_pd( SSE2_MUL_CONJ( *U1 , *S1 ) ,
				  SSE2_MULCONJ( *U2 , *S2 ) ) ) ;
    sm1 = _mm_add_pd( sm1 ,
		      _mm_sub_pd( SSE2_MUL_CONJ( *U1 , *S2 ) ,
				  SSE2_MULCONJ( *U2 , *S1 ) ) ) ;
    // increment our pointers
    S1++ , S2++ , U1++ , U2++ ;
  }
#endif

  // puts the norm in both parts
  register __m128d z = _mm_add_pd( _mm_mul_pd( sm0 , sm0 ) ,
				   _mm_mul_pd( sm1 , sm1 ) ) ;
  z = _mm_add_pd( z , _mm_shuffle_pd( z , z , 1 ) ) ;
  z = _mm_sqrt_pd( z ) ;
  z = _mm_div_pd( _mm_set1_pd( 1.0 ) , z ) ;
  sm0 = _mm_mul_pd( sm0 , z ) ;
  sm1 = _mm_mul_pd( sm1 , z ) ;

  // poke back into *s0 and *s1 and *scale
  _mm_store_pd( (void*)s0 , sm0 ) ; 
  _mm_store_pd( (void*)s1 , sm1 ) ; 
  _mm_store_sd( (void*)scale , z ) ;

  return ;
}

// compactified (sparse matrix rep) su(2) multiply of,
//
//     | a b || M[row(a)] M[row(a)++] .... |   
//     | c d || M[row(c)] M[row(c)++] .... |
//
void
shortened_su2_multiply( GLU_complex *w , 
			const GLU_complex a , 
			const GLU_complex b , 
			const GLU_complex c , 
			const GLU_complex d , 
			const size_t su2_index )
{
  const register __m128d A = _mm_setr_pd( creal( a ) , cimag( a ) ) ;
  const register __m128d B = _mm_setr_pd( creal( b ) , cimag( b ) ) ;

  const size_t row_a = NC * (int)( Latt.su2_data[ su2_index ].idx_a / NC ) ;
  const size_t row_c = NC * (int)( Latt.su2_data[ su2_index ].idx_c / NC ) ;

  register __m128d tmp ;
  __m128d *w1 = (__m128d*)( w + row_a ) ;
  __m128d *w2 = (__m128d*)( w + row_c ) ;

  size_t i ;
  for( i = 0 ; i < NC ; i++ ) {
    tmp = *w1 ;
    *w1 = _mm_add_pd( SSE2_MUL( A , *w1 ) , SSE2_MUL( B , *w2 ) ) ;
    *w2 = _mm_sub_pd( SSE2_MULCONJ( A , *w2 ) , SSE2_MULCONJ( B , tmp ) ) ;
    w1++ , w2++ ;
  }
  return ;
}

//  compactified M.su(2)^{\dagger} multiply of,
//
//   | M[col(a)]    M[col(b)]    | | a b |^{\dagger}
//   | M[col(a)+NC] M[col(b)+NC] | | c d |
//   | .....        .......      |
//
void
shortened_su2_multiply_dag( GLU_complex *U , 
			    const GLU_complex a , 
			    const GLU_complex b , 
			    const GLU_complex c , 
			    const GLU_complex d , 
			    const size_t su2_index )
{
  GLU_complex U1 , U2 ; // temporaries for caching
  const size_t col_a = (int)( Latt.su2_data[ su2_index ].idx_a % NC ) ;
  const size_t col_b = (int)( Latt.su2_data[ su2_index ].idx_b % NC ) ;
  size_t i ;
  for( i = 0 ; i < NC ; i++ ) {
    U1 = U[ col_a + i*NC ] ;
    U2 = U[ col_b + i*NC ] ;
    U[ col_a + i*NC ] = U1 * conj(a) + U2 * conj(b) ;
    U[ col_b + i*NC ] = U1 * conj(c) + U2 * conj(d) ;
  }
  return ;
}

// rotate a matrix U = su2_i*U where su2_i is an su2 matrix embedded in suN
void
su2_rotate( GLU_complex U[ NCNC ] ,
	    const GLU_complex s0 ,
	    const GLU_complex s1 ,
	    const size_t su2_index )
{
#if NC == 3
  __m128d *u = (__m128d*)U ;
  register const __m128d sm0 = _mm_setr_pd( creal( s0 ) , cimag( s0 ) ) ;
  register const __m128d sm1 = _mm_setr_pd( creal( s1 ) , cimag( s1 ) ) ;
  register __m128d tmp0 , tmp1 , a , b ;
  switch( su2_index%3 ) { // again I don't like this
  case 0 :
    // first one
    a = *( u + 0 ) ; b = *( u + 3 ) ;
    tmp0 = _mm_add_pd( SSE2_MUL( sm0 , a ) ,
		       SSE2_MUL( sm1 , b ) ) ;
    tmp1 = _mm_sub_pd( SSE2_MULCONJ( sm0 , b ) ,
		       SSE2_MULCONJ( sm1 , a ) ) ;
    *( u + 0 ) = tmp0 ;
    *( u + 3 ) = tmp1 ;
    // second one
    a = *( u + 1 ) ; b = *( u + 4 ) ;
    tmp0 = _mm_add_pd( SSE2_MUL( sm0 , a ) ,
		       SSE2_MUL( sm1 , b ) ) ;
    tmp1 = _mm_sub_pd( SSE2_MULCONJ( sm0 , b ) ,
		       SSE2_MULCONJ( sm1 , a ) ) ;
    *( u + 1 ) = tmp0 ;
    *( u + 4 ) = tmp1 ;
    // third
    a = *( u + 2 ) ; b = *( u + 5 ) ;
    tmp0 = _mm_add_pd( SSE2_MUL( sm0 , a ) ,
		       SSE2_MUL( sm1 , b ) ) ;
    tmp1 = _mm_sub_pd( SSE2_MULCONJ( sm0 , b ) ,
		       SSE2_MULCONJ( sm1 , a ) ) ;
    *( u + 2 ) = tmp0 ;
    *( u + 5 ) = tmp1 ;
    break ;
  case 1 :
    // first one
    a = *( u + 3 ) ; b = *( u + 6 ) ;
    tmp0 = _mm_add_pd( SSE2_MUL( sm0 , a ) ,
		       SSE2_MUL( sm1 , b ) ) ;
    tmp1 = _mm_sub_pd( SSE2_MULCONJ( sm0 , b ) ,
		       SSE2_MULCONJ( sm1 , a ) ) ;
    *( u + 3 ) = tmp0 ;
    *( u + 6 ) = tmp1 ;
    // second one
    a = *( u + 4 ) ; b = *( u + 7 ) ;
    tmp0 = _mm_add_pd( SSE2_MUL( sm0 , a ) ,
		       SSE2_MUL( sm1 , b ) ) ;
    tmp1 = _mm_sub_pd( SSE2_MULCONJ( sm0 , b ) ,
		       SSE2_MULCONJ( sm1 , a ) ) ;
    *( u + 4 ) = tmp0 ;
    *( u + 7 ) = tmp1 ;
    // third
    a = *( u + 5 ) ; b = *( u + 8 ) ;
    tmp0 = _mm_add_pd( SSE2_MUL( sm0 , a ) ,
		       SSE2_MUL( sm1 , b ) ) ;
    tmp1 = _mm_sub_pd( SSE2_MULCONJ( sm0 , b ) ,
		       SSE2_MULCONJ( sm1 , a ) ) ;
    *( u + 5 ) = tmp0 ;
    *( u + 8 ) = tmp1 ;
    break ;
  case 2 :
    // first one
    a = *( u + 0 ) ; b = *( u + 6 ) ;
    tmp0 = _mm_sub_pd( SSE2_MULCONJ( sm0 , a ) ,
		       SSE2_MULCONJ( sm1 , b ) ) ;
    tmp1 = _mm_add_pd( SSE2_MUL( sm0 , b ) ,
		       SSE2_MUL( sm1 , a ) ) ;
    *( u + 0 ) = tmp0 ;
    *( u + 6 ) = tmp1 ;
    // second
    a = *( u + 1 ) ; b = *( u + 7 ) ;
    tmp0 = _mm_sub_pd( SSE2_MULCONJ( sm0 , a ) ,
		       SSE2_MULCONJ( sm1 , b ) ) ;
    tmp1 = _mm_add_pd( SSE2_MUL( sm0 , b ) ,
		       SSE2_MUL( sm1 , a ) ) ;
    *( u + 1 ) = tmp0 ;
    *( u + 7 ) = tmp1 ;
    // third
    a = *( u + 2 ) ; b = *( u + 8 ) ;
    tmp0 = _mm_sub_pd( SSE2_MULCONJ( sm0 , a ) ,
		       SSE2_MULCONJ( sm1 , b ) ) ;
    tmp1 = _mm_add_pd( SSE2_MUL( sm0 , b ) ,
		       SSE2_MUL( sm1 , a ) ) ;
    *( u + 2 ) = tmp0 ;
    *( u + 8 ) = tmp1 ;
    break ;
  }
#elif NC == 2
  __m128d *u = (__m128d*)U ;
  register const __m128d sm0 = _mm_setr_pd( creal( s0 ) , cimag( s0 ) ) ;
  register const __m128d sm1 = _mm_setr_pd( creal( s1 ) , cimag( s1 ) ) ;
  *( u + 0 ) = _mm_add_pd( SSE2_MUL( sm0 , *( u + 0 ) ) ,
			   SSE2_MUL( sm1 , *( u + 2 ) ) ) ;
  *( u + 1 ) = _mm_add_pd( SSE2_MUL( sm0 , *( u + 1 ) ) ,
			   SSE2_MUL( sm1 , *( u + 3 ) ) ) ;
  *( u + 2 ) = SSE_FLIP( SSE2_CONJ( *( u + 1 ) ) ) ; 
  *( u + 3 ) = SSE2_CONJ( *( u + 0 ) ) ;
#else
  // just a call to su2 multiply
  shortened_su2_multiply( U , s0 , s1 , -conj(s1) , conj(s0) , su2_index ) ;
#endif
  return ;
}

#endif
