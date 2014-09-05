/*
    Copyright 2013 Renwick James Hudspith

    This file (givens.c) is part of GLU.

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
   @file givens.c
   @brief method based on givens rotation for trace maximisation

   Some credit to Urs Wenger and ETMC and Edwards of Chroma fame

   performs the multiplication
   w -> su(2)_i * w
   U -> U * su(2)_i 

   @warning single precision code with gcc-4.7.2 cannot compile this, I have no idea why: is known bug in gcc
 */

#include "Mainfile.h"

/**
   @param GIVE_PRECOND
   @brief precondition the routine by giving the reunitarised matrix as an initial guess
   This is seen to reduce the number of iterations and preserve gauge covariance in the APE and HYP smearing for large #NC
 */
#define GIVE_PRECOND

#ifdef GIVE_PRECOND
  #include "gramschmidt.h" 
#endif

// allocate the struct somewhere
static struct su2_subgroups *su2_data ;

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
			const int su2_index )
{
  GLU_complex W1 , W2 ; // temporaries
  const int row_a = NC * (int)( su2_data[ su2_index ].idx_a / NC ) ;
  const int row_c = NC * (int)( su2_data[ su2_index ].idx_c / NC ) ;
  int i ;
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
			    const int su2_index )
{
  GLU_complex U1 , U2 ; // temporaries for caching
  const int col_a = (int)( su2_data[ su2_index ].idx_a % NC ) ;
  const int col_b = (int)( su2_data[ su2_index ].idx_b % NC ) ;
  int i ;
  for( i = 0 ; i < NC ; i++ ) {
    U1 = U[ col_a + i*NC ] ;
    U2 = U[ col_b + i*NC ] ;
    U[ col_a + i*NC ] = U1 * conj(a) + U2 * conj(b) ;
    U[ col_b + i*NC ] = U1 * conj(c) + U2 * conj(d) ;
  }
  return ;
}

// overrelaxation
static void
overrelax( GLU_complex *s0 , 
	   GLU_complex *s1 , 
	   const double OrParam )
{
  // I use the MILC OR step here
  register const GLU_real asq = cimag(*s0)*cimag(*s0) 
    + creal(*s1)*creal(*s1) + cimag(*s1)*cimag(*s1) ;
  register const GLU_real a0sq = creal(*s0)*creal(*s0) ;
  register const GLU_real x = ( OrParam * a0sq + asq ) / ( a0sq + asq ) ;
  register const GLU_real r = sqrt( a0sq + x*x*asq ) ;
  register const GLU_real xdr = x/r ;
  *s0 = creal( *s0 ) / r + I * cimag( *s0 ) * xdr ;
  *s1 *= xdr ;
  return ;
}

#if NC == 3

//  The first of the su(2) subgroups for su(3) is
//
//  |  s0   s1  0 |
//  | -s1*  s0* 0 |
//  |  0     0  1 |
//
static inline void
rotation1( U , w )
     GLU_complex U[ NCNC ] , w[ NCNC ] ;
{
  register GLU_complex s0 = U[0] + conj( U[4] ) ;
  register GLU_complex s1 = U[1] - conj( U[3] ) ;
  const double scale = 1.0 / sqrt( creal(s0)*creal(s0) + cimag(s0)*cimag(s0) + \
				   creal(s1)*creal(s1) + cimag(s1)*cimag(s1) ) ;
  s0 *= scale ;
  s1 *= scale ;

  // complex conjugates
  register GLU_complex s0_star = conj( s0 ) ;
  register GLU_complex s1_star = conj( s1 ) ;
  // compute the overwritten matrix "w -> su(2)_1 * w"
  GLU_complex tmp0 = s0 * w[0] + s1 * w[3] ;
  GLU_complex tmp1 = s0 * w[1] + s1 * w[4] ;
  GLU_complex tmp2 = s0 * w[2] + s1 * w[5] ;
  GLU_complex tmp3 = -s1_star * w[0] + s0_star * w[3] ;
  GLU_complex tmp4 = -s1_star * w[1] + s0_star * w[4] ;
  GLU_complex tmp5 = -s1_star * w[2] + s0_star * w[5] ;
  *( w + 0 ) = tmp0 ;
  *( w + 1 ) = tmp1 ;
  *( w + 2 ) = tmp2 ;
  *( w + 3 ) = tmp3 ;
  *( w + 4 ) = tmp4 ;
  *( w + 5 ) = tmp5 ;
  // compute the new matrix "U -> U * su(2)_1" these are columnal!
  tmp0 = U[0] * s0_star + U[1] * s1_star ;
  tmp1 = -U[0] * s1 + U[1] * s0 ;
  tmp2 = U[3] * s0_star + U[4] * s1_star ;
  tmp3 = -U[3] * s1 + U[4] * s0 ;
  tmp4 = U[6] * s0_star + U[7] * s1_star ;
  tmp5 = -U[6] * s1 + U[7] * s0 ;
  *( U + 0 ) = tmp0 ;
  *( U + 1 ) = tmp1 ;
  *( U + 3 ) = tmp2 ;
  *( U + 4 ) = tmp3 ;
  *( U + 6 ) = tmp4 ;
  *( U + 7 ) = tmp5 ;
  return ;
}

// second rotation matrix is
//
//  |  1    0   0  |
//  |  0   s0  s1  |
//  |  0  -s1* s0* |
//
static inline void
rotation2( U , w )
     GLU_complex U[ NCNC ] , w[ NCNC ] ;
{
  register GLU_complex s0 = U[4] + conj( U[8] ) ;
  register GLU_complex s1 = U[5] - conj( U[7] ) ;
  const double scale = 1.0 / sqrt( creal(s0)*creal(s0) + cimag(s0)*cimag(s0) + \
				   creal(s1)*creal(s1) + cimag(s1)*cimag(s1) ) ;
  s0 *= scale ;
  s1 *= scale ;

  // complex conjugates
  register const GLU_complex s0_star = conj( s0 ) ;
  register const GLU_complex s1_star = conj( s1 ) ;
  // compute the overwritten matrix "w -> su(2)_1 * w"
  GLU_complex tmp0 = s0 * w[3] + s1 * w[6] ;
  GLU_complex tmp1 = s0 * w[4] + s1 * w[7] ;
  GLU_complex tmp2 = s0 * w[5] + s1 * w[8] ;
  GLU_complex tmp3 = -s1_star * w[3] + s0_star * w[6] ;
  GLU_complex tmp4 = -s1_star * w[4] + s0_star * w[7] ;
  GLU_complex tmp5 = -s1_star * w[5] + s0_star * w[8] ;
  *( w + 3 ) = tmp0 ;
  *( w + 4 ) = tmp1 ;
  *( w + 5 ) = tmp2 ;
  *( w + 6 ) = tmp3 ;
  *( w + 7 ) = tmp4 ;
  *( w + 8 ) = tmp5 ;
  // compute the new matrix "U -> U * su(2)_1" these are columnal!
  tmp0 = U[1] * s0_star + U[2] * s1_star ;
  tmp1 = -U[1] * s1 + U[2] * s0 ;
  tmp2 = U[4] * s0_star + U[5] * s1_star ;
  tmp3 = -U[4] * s1 + U[5] * s0 ;
  tmp4 = U[7] * s0_star + U[8] * s1_star ;
  tmp5 = -U[7] * s1 + U[8] * s0 ;
  *( U + 1 ) = tmp0 ;
  *( U + 2 ) = tmp1 ;
  *( U + 4 ) = tmp2 ;
  *( U + 5 ) = tmp3 ;
  *( U + 7 ) = tmp4 ;
  *( U + 8 ) = tmp5 ;
  return ;
}

// The third rotation looks like
//
//  | s0  0  s1   |
//  |  0  1  0    |
//  | -s1* 0  s0* |
//
static inline void
rotation3( U , w )
     GLU_complex U[ NCNC ] , w[ NCNC ] ;
{
  register GLU_complex s0 = U[8] + conj( U[0] ) ;
  register GLU_complex s1 = U[6] - conj( U[2] ) ;
  const double scale = 1.0 / sqrt( creal(s0)*creal(s0) + cimag(s0)*cimag(s0) + \
				   creal(s1)*creal(s1) + cimag(s1)*cimag(s1) ) ;
  s0 *= scale ;
  s1 *= scale ;
  // complex conjugates
  register const GLU_complex s0_star = conj( s0 ) ;
  register const GLU_complex s1_star = conj( s1 ) ;
  // compute the overwritten matrix "w -> su(2)_1 * w"
  GLU_complex tmp0 = s0_star * w[0] - s1_star * w[6] ;
  GLU_complex tmp1 = s0_star * w[1] - s1_star * w[7] ;
  GLU_complex tmp2 = s0_star * w[2] - s1_star * w[8] ;
  GLU_complex tmp3 = s1 * w[0] + s0 * w[6] ;
  GLU_complex tmp4 = s1 * w[1] + s0 * w[7] ;
  GLU_complex tmp5 = s1 * w[2] + s0 * w[8] ;
  *( w + 0 ) = tmp0 ;
  *( w + 1 ) = tmp1 ;
  *( w + 2 ) = tmp2 ;
  *( w + 6 ) = tmp3 ;
  *( w + 7 ) = tmp4 ;
  *( w + 8 ) = tmp5 ;
  tmp0 = U[0] * s0 - U[2] * s1 ;
  tmp1 = U[0] * s1_star + U[2] * s0_star ;
  tmp2 = U[3] * s0 - U[5] * s1 ;
  tmp3 = U[3] * s1_star + U[5] * s0_star ;
  tmp4 = U[6] * s0 - U[8] * s1 ;
  tmp5 = U[6] * s1_star + U[8] * s0_star ;
  *( U + 0 ) = tmp0 ;
  *( U + 2 ) = tmp1 ;
  *( U + 3 ) = tmp2 ;
  *( U + 5 ) = tmp3 ;
  *( U + 6 ) = tmp4 ;
  *( U + 8 ) = tmp5 ;
  return ;
}

#else

// NC generic givens rotations
static void
rotation( U , w , su2_index )
     GLU_complex U[ NCNC ] , w[ NCNC ] ;
     const int su2_index ;
{
  const int a = su2_data[ su2_index ].idx_a ;
  const int b = su2_data[ su2_index ].idx_b ;
  const int c = su2_data[ su2_index ].idx_c ;
  const int d = su2_data[ su2_index ].idx_d ;

  register GLU_complex s0 = U[a] + conj( U[d] ) ;
  register GLU_complex s1 = U[b] - conj( U[c] ) ;

  // I know how to hypot this, but don't know whether it is worth it
  const double scale = 1.0 / sqrt( creal(s0)*creal(s0) + cimag(s0)*cimag(s0) + \
				   creal(s1)*creal(s1) + cimag(s1)*cimag(s1) ) ;
  s0 *= scale ;
  s1 *= scale ;

  // these overwrite w and U, utilising the fact that most 
  // of the rotation matrixis filled with zeros, ONLY does 
  // the necessary operations.
  // Factor of 10x for SU(8) for doing this, old code is back 
  // in the SVN for checking
  shortened_su2_multiply( w , s0 , s1 , -conj(s1) , conj(s0) , su2_index ) ;
  shortened_su2_multiply_dag( U , s0 , s1 , -conj(s1) , conj(s0) , su2_index ) ;

  return ;
}
#endif

// compute the relevant su2 indices
void
compute_pertinent_indices( void )
{
  su2_data = (struct su2_subgroups*)malloc( NSU2SUBGROUPS * sizeof( struct su2_subgroups ) ) ;
  int su2_index ;
  for( su2_index = 0 ; su2_index < NSU2SUBGROUPS ; su2_index ++ ) {
    int i1, i2 ;
    int found = 0 ;
    int del_i = 0 ;
    int index = -1 ;
    
    while ( del_i < (NC-1) && found == 0 ) {
      del_i++;
      for ( i1 = 0; i1 < (NC-del_i); i1++ ) {
	index++;
	if ( index == su2_index ) {
	  found = 1;
	  break;
	}
      }
    }
    i2 = i1 + del_i;
    
    su2_data[ su2_index ].idx_a = i1 + NC * i1 ;
    su2_data[ su2_index ].idx_b = i2 + NC * i1 ;
    su2_data[ su2_index ].idx_c = i1 + NC * i2 ;
    su2_data[ su2_index ].idx_d = i2 + NC * i2 ;
  }
  return ;
}

// free the struct
void
free_su2_data()
{
  free( su2_data ) ;
}

// Cabibbo-Marinari projections for trace maximisation
void
givens_reunit( GLU_complex U[ NCNC ] ) 
{
  GLU_complex w[ NCNC ] ;
  double trace_new , trace_old ;

  #ifdef GIVE_PRECOND
  // MILC precondition,
  // provides a starting guess; the reunitarised path. This is seen 
  // to consistently reduce the number of iterations needed.
  // Seems to work (i.e. reduce breaking of gauge invariance
  // for APE smearing for large NC) Big slow down though due to gram 
  // schmidt and nested SU(2) updates; for SU(3) looks pretty good though 
  // suggest using n-ape for large NC though, is much stabler
  GLU_complex temp[ NCNC ] ;
  equiv( temp , U ) ;
  equiv( w , U ) ;
  reunit2( w ) ;
  multab_dag( U , temp , w ) ;
  #else
  identity( w ) ;
  #endif

  speed_trace_Re( &trace_old , U ) ;

  int i ;
  // seventy five iterations is quite high no? Code overwrites
  // U with the product U.W^\dagger
  for( i = 0 ; i < 75 ; i++ ) {

    #if NC == 3

    rotation1( U , w ) ;
    rotation2( U , w ) ;
    rotation3( U , w ) ;

    trace_new = (double)creal( U[0] ) 
              + (double)creal( U[4] ) 
              + (double)creal( U[8] ) ;

    #elif NC == 2

    GLU_complex rot[ NCNC ] , tmp[ NCNC ] ;
    register GLU_complex s0 = U[0] + conj( U[3] ) ;
    register GLU_complex s1 = U[1] - conj( U[2] ) ;
    const double scale = 1.0 / sqrt( creal(s0)*creal(s0) + cimag(s0)*cimag(s0) + \
				     creal(s1)*creal(s1) + cimag(s1)*cimag(s1) ) ;
    s0 *= scale ;
    s1 *= scale ;
    rot[0] = s0 ;
    rot[3] = conj( s0 ) ;
    rot[1] = s1 ;
    rot[2] = -conj( s1 ) ;
    multab_suNC( tmp , rot , w ) ; // is guaranteed SU(N)
    equiv( w , tmp ) ;
    multab_dag( tmp , U , rot ) ; // is not
    equiv( U , tmp ) ;
    trace_new = (double)creal( U[0] ) + (double)creal( U[3] ) ;

    // su2 only ever needs one update ...                                       
    break ;

    #else

    // NC-generic cabbibo-marinari update over the su2 subgroups
    int su2_index ;
    for( su2_index = 0 ; su2_index < NSU2SUBGROUPS ; su2_index++ ) {
      rotation( U , w , su2_index ) ;
    }      
    speed_trace_Re( &trace_new , U ) ;
 
    #endif
    
    // check if the newest trace product is maximised
    if( fabs( trace_new - trace_old ) / ( trace_new ) < PREC_TOL ) { break ; }
    trace_old = trace_new ;
  }
  // set U == to the maximal trace-inducing rotation
  equiv( U , w ) ;
  return ;
}

// NC generic givens rotations
void
OrRotation( const GLU_complex U[ NCNC ] , 
	    GLU_complex *s0 , 
	    GLU_complex *s1 ,
	    const double OrParam ,
	    const int su2_index )
{
  *s0 = U[su2_data[ su2_index ].idx_a] 
    + conj( U[su2_data[ su2_index ].idx_d] ) ;
  *s1 = U[su2_data[ su2_index ].idx_b] 
    - conj( U[su2_data[ su2_index ].idx_c] ) ;

  // overrelax
  overrelax( s0 , s1 , OrParam ) ;
  return ;
}

// clean up for local scope
#ifdef GIVE_PRECOND
  #undef GIVE_PRECOND
#endif

