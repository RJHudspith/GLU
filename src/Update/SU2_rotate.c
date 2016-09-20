/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (SU2_rotate.c) is part of GLU.

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
   @file SU2_rotate.c
   @brief SU2 subgroup rotations on SU(#NC) matrices
 */
#include "Mainfile.h"

#if !( defined HAVE_IMMINTRIN_H ) || ( defined SINGLE_PREC )

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
  register GLU_complex temp0 = 0.0 , temp1 = 0.0 , temp2 = 0.0 , temp3 = 0.0 ;
#if NC == 3
  switch( su2_index%3 ) {
    // rotation 1
    //  |  s0   s1  0 |
    //  | -s1*  s0* 0 |
    //  |  0     0  1 |
  case 0 :
    temp0 = U[0]*conj(staple[0])+U[1]*conj(staple[1])+U[2]*conj(staple[2]) ;
    temp1 = U[0]*conj(staple[3])+U[1]*conj(staple[4])+U[2]*conj(staple[5]) ;
    temp2 = U[3]*conj(staple[0])+U[4]*conj(staple[1])+U[5]*conj(staple[2]) ;
    temp3 = U[3]*conj(staple[3])+U[4]*conj(staple[4])+U[5]*conj(staple[5]) ;
    break ;
  case 1 :
    // rotation 2
    //  |  1    0   0  |
    //  |  0   s0  s1  |
    //  |  0  -s1* s0* |
    temp0 = U[3]*conj(staple[3])+U[4]*conj(staple[4])+U[5]*conj(staple[5]) ;
    temp1 = U[3]*conj(staple[6])+U[4]*conj(staple[7])+U[5]*conj(staple[8]) ;
    temp2 = U[6]*conj(staple[3])+U[7]*conj(staple[4])+U[8]*conj(staple[5]) ;
    temp3 = U[6]*conj(staple[6])+U[7]*conj(staple[7])+U[8]*conj(staple[8]) ;
    break ;
  case 2 :
    // rotation 3
    //  | s0*  0  -s1 |
    //  |  0   1   0  |
    //  | s1   0  s0  |
    temp0 = U[6]*conj(staple[6])+U[7]*conj(staple[7])+U[8]*conj(staple[8]) ;
    temp1 = U[6]*conj(staple[0])+U[7]*conj(staple[1])+U[8]*conj(staple[2]) ;
    temp2 = U[0]*conj(staple[6])+U[1]*conj(staple[7])+U[2]*conj(staple[8]) ;
    temp3 = U[0]*conj(staple[0])+U[1]*conj(staple[1])+U[2]*conj(staple[2]) ;
    break ;
  }
#elif NC == 2
  temp0 = U[0] * conj( staple[0] ) + U[1] * conj( staple[1] ) ;
  temp1 = U[0] * conj( staple[2] ) + U[1] * conj( staple[3] ) ;
  temp2 = U[2] * conj( staple[0] ) + U[3] * conj( staple[1] ) ;
  temp3 = U[2] * conj( staple[2] ) + U[3] * conj( staple[3] ) ;
#else
  // su(N) version
  const size_t row_a = Latt.su2_data[ su2_index ].idx_a / NC ;
  const size_t col_b = Latt.su2_data[ su2_index ].idx_b % NC ;
  size_t i ;
  for( i = 0 ; i < NC ; i++ ) {
    temp0 += U[ i + NC * row_a ] * conj( staple[ i + NC * row_a ] ) ;
    temp1 += U[ i + NC * row_a ] * conj( staple[ i + NC * col_b ] ) ;
    temp2 += U[ i + NC * col_b ] * conj( staple[ i + NC * row_a ] ) ;
    temp3 += U[ i + NC * col_b ] * conj( staple[ i + NC * col_b ] ) ;
  }
#endif
  *s0 = temp0 + conj( temp3 ) ;
  *s1 = temp1 - conj( temp2 ) ;
  // pass the inverse determinant
  *scale =								\
    1.0 / sqrt( creal(*s0)*creal(*s0) + cimag(*s0)*cimag(*s0) +		\
		creal(*s1)*creal(*s1) + cimag(*s1)*cimag(*s1) ) ;
  *s0 *= (*scale) ;
  *s1 *= (*scale) ;
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
  GLU_complex W1 , W2 ; // temporaries
  const size_t row_a = NC * (int)( Latt.su2_data[ su2_index ].idx_a / NC ) ;
  const size_t row_c = NC * (int)( Latt.su2_data[ su2_index ].idx_c / NC ) ;
  size_t i ;
  for( i = 0 ; i < NC ; i++ ) {
    W1 = w[ row_a + i ] ;
    W2 = w[ row_c + i ] ;
    w[ row_a + i ] = a * W1 + b * W2 ;
    w[ row_c + i ] = c * W1 + d * W2 ;
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
  register GLU_complex s0_star = conj( s0 ) ;
  register GLU_complex s1_star = conj( s1 ) ;
  register GLU_complex tmp0 , tmp1 , tmp2 , tmp3 , tmp4 , tmp5 ;
  switch( su2_index%3 ) {
  case 0 :
    tmp0 = s0 * U[0] + s1 * U[3] ;
    tmp1 = s0 * U[1] + s1 * U[4] ;
    tmp2 = s0 * U[2] + s1 * U[5] ;
    tmp3 = -s1_star * U[0] + s0_star * U[3] ;
    tmp4 = -s1_star * U[1] + s0_star * U[4] ;
    tmp5 = -s1_star * U[2] + s0_star * U[5] ;
    *( U + 0 ) = tmp0 ;
    *( U + 1 ) = tmp1 ;
    *( U + 2 ) = tmp2 ;
    *( U + 3 ) = tmp3 ;
    *( U + 4 ) = tmp4 ;
    *( U + 5 ) = tmp5 ;
    break ;
  case 1 :
    tmp0 = s0 * U[3] + s1 * U[6] ;
    tmp1 = s0 * U[4] + s1 * U[7] ;
    tmp2 = s0 * U[5] + s1 * U[8] ;
    tmp3 = -s1_star * U[3] + s0_star * U[6] ;
    tmp4 = -s1_star * U[4] + s0_star * U[7] ;
    tmp5 = -s1_star * U[5] + s0_star * U[8] ;
    *( U + 3 ) = tmp0 ;
    *( U + 4 ) = tmp1 ;
    *( U + 5 ) = tmp2 ;
    *( U + 6 ) = tmp3 ;
    *( U + 7 ) = tmp4 ;
    *( U + 8 ) = tmp5 ;
    break ;
  case 2 :
    tmp0 = s0_star * U[0] - s1_star * U[6] ;
    tmp1 = s0_star * U[1] - s1_star * U[7] ;
    tmp2 = s0_star * U[2] - s1_star * U[8] ;
    tmp3 = s1 * U[0] + s0 * U[6] ;
    tmp4 = s1 * U[1] + s0 * U[7] ;
    tmp5 = s1 * U[2] + s0 * U[8] ;
    *( U + 0 ) = tmp0 ;
    *( U + 1 ) = tmp1 ;
    *( U + 2 ) = tmp2 ;
    *( U + 6 ) = tmp3 ;
    *( U + 7 ) = tmp4 ;
    *( U + 8 ) = tmp5 ;
    break ;
  }
#elif NC == 2
  // U = s * U
  *( U + 0 ) = s0 * U[0] + s1 * U[2] ;
  *( U + 1 ) = s0 * U[1] + s1 * U[3] ;
  *( U + 2 ) = -conj( U[1] ) ;
  *( U + 3 ) =  conj( U[0] ) ;
#else
  // just a call to su2 multiply
  shortened_su2_multiply( U , s0 , s1 , -conj(s1) , conj(s0) , su2_index ) ;
#endif
  return ;
}

#endif
