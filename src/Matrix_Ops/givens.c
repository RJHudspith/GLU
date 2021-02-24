/*
    Copyright 2013-2018 Renwick James Hudspith

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

#include "givens.h"     // alphabetising
#include "SU2_rotate.h" // su2 subgroup rotations

/**
   @param GIVE_PRECOND
   @brief precondition the routine by giving the reunitarised matrix as an initial guess
   This is seen to reduce the number of iterations and preserve gauge covariance in the APE and HYP smearing for large #NC
 */
#define GIVE_PRECOND

#ifdef GIVE_PRECOND
  #include "gramschmidt.h" 
#endif

#if NC == 3

#if (defined HAVE_IMMINTRIN_H) && !(defined SINGLE_PREC)

#include <immintrin.h>
#include "SSE2_OPS.h"

static inline void
rotation1( GLU_complex *U , 
	   GLU_complex *w )
{
  __m128d *pU = (__m128d*)U ;
  __m128d *pw = (__m128d*)w ;
  register __m128d s0 = _mm_add_pd( *(pU + 0) , SSE2_CONJ( *(pU + 4) ) ) ;
  register __m128d s1 = _mm_sub_pd( *(pU + 1) , SSE2_CONJ( *(pU + 3) ) ) ;
  register __m128d sum = _mm_add_pd( _mm_mul_pd( s0 , s0 ) ,
				     _mm_mul_pd( s1 , s1 ) ) ;
  sum = _mm_add_pd( sum , _mm_shuffle_pd( sum , sum , 1 ) ) ;
  sum = _mm_sqrt_pd( sum ) ;
  s0  = _mm_div_pd( s0 , sum ) ;
  s1  = _mm_div_pd( s1 , sum ) ;  
  // compute the overwritten matrix "w -> su(2)_1 * w"
  __m128d tmp0 = _mm_add_pd( SSE2_MUL( s0 , *(pw + 0) ) ,
			     SSE2_MUL( s1 , *(pw + 3) ) ) ;
  __m128d tmp1 = _mm_add_pd( SSE2_MUL( s0 , *(pw + 1) ) ,
			     SSE2_MUL( s1 , *(pw + 4) ) ) ;
  __m128d tmp2 = _mm_add_pd( SSE2_MUL( s0 , *(pw + 2) ) ,
			     SSE2_MUL( s1 , *(pw + 5) ) ) ;
  __m128d tmp3 = _mm_sub_pd( SSE2_MULCONJ( s0 , *(pw + 3) ) ,
			     SSE2_MULCONJ( s1 , *(pw + 0) ) ) ;
  __m128d tmp4 = _mm_sub_pd( SSE2_MULCONJ( s0 , *(pw + 4) ) ,
			     SSE2_MULCONJ( s1 , *(pw + 1) ) ) ;
  __m128d tmp5 = _mm_sub_pd( SSE2_MULCONJ( s0 , *(pw + 5) ) ,
			     SSE2_MULCONJ( s1 , *(pw + 2) ) ) ;	     
  *( pw + 0 ) = tmp0 ;
  *( pw + 1 ) = tmp1 ;
  *( pw + 2 ) = tmp2 ;
  *( pw + 3 ) = tmp3 ;
  *( pw + 4 ) = tmp4 ;
  *( pw + 5 ) = tmp5 ;
  // compute the new matrix "U -> U * su(2)_1" these are columnal!
  tmp0 = _mm_add_pd( SSE2_MULCONJ( s0 , *(pU + 0) ) ,
		     SSE2_MULCONJ( s1 , *(pU + 1) ) ) ;
  tmp1 = _mm_sub_pd( SSE2_MUL( s0 , *(pU + 1) ) ,
		     SSE2_MUL( s1 , *(pU + 0) ) ) ;
  tmp2 = _mm_add_pd( SSE2_MULCONJ( s0 , *(pU + 3) ) ,
		     SSE2_MULCONJ( s1 , *(pU + 4) ) ) ;
  tmp3 = _mm_sub_pd( SSE2_MUL( s0 , *(pU + 4) ) ,
		     SSE2_MUL( s1 , *(pU + 3) ) ) ;
  tmp4 = _mm_add_pd( SSE2_MULCONJ( s0 , *(pU + 6) ) ,
		     SSE2_MULCONJ( s1 , *(pU + 7) ) ) ;
  tmp5 = _mm_sub_pd( SSE2_MUL( s0 , *(pU + 7) ) ,
		     SSE2_MUL( s1 , *(pU + 6) ) ) ;
  *( pU + 0 ) = tmp0 ;
  *( pU + 1 ) = tmp1 ;
  *( pU + 3 ) = tmp2 ;
  *( pU + 4 ) = tmp3 ;
  *( pU + 6 ) = tmp4 ;
  *( pU + 7 ) = tmp5 ;
  return ;
}

static inline void
rotation2( GLU_complex *U , 
	   GLU_complex *w )
{
  __m128d *pU = (__m128d*)U ;
  __m128d *pw = (__m128d*)w ;
  register __m128d s0 = _mm_add_pd( *(pU + 4) , SSE2_CONJ( *(pU + 8) ) ) ;
  register __m128d s1 = _mm_sub_pd( *(pU + 5) , SSE2_CONJ( *(pU + 7) ) ) ;
  register __m128d sum = _mm_add_pd( _mm_mul_pd( s0 , s0 ) ,
				     _mm_mul_pd( s1 , s1 ) ) ;
  sum = _mm_add_pd( sum , _mm_shuffle_pd( sum , sum , 1 ) ) ;
  sum = _mm_sqrt_pd( sum ) ;
  s0  = _mm_div_pd( s0 , sum ) ;
  s1  = _mm_div_pd( s1 , sum ) ;  
  // compute the overwritten matrix "w -> su(2)_1 * w"
  __m128d tmp0 = _mm_add_pd( SSE2_MUL( s0 , *(pw + 3) ) ,
			     SSE2_MUL( s1 , *(pw + 6) ) ) ;
  __m128d tmp1 = _mm_add_pd( SSE2_MUL( s0 , *(pw + 4) ) ,
			     SSE2_MUL( s1 , *(pw + 7) ) ) ;
  __m128d tmp2 = _mm_add_pd( SSE2_MUL( s0 , *(pw + 5) ) ,
			     SSE2_MUL( s1 , *(pw + 8) ) ) ;
  __m128d tmp3 = _mm_sub_pd( SSE2_MULCONJ( s0 , *(pw + 6) ) ,
			     SSE2_MULCONJ( s1 , *(pw + 3) ) ) ;
  __m128d tmp4 = _mm_sub_pd( SSE2_MULCONJ( s0 , *(pw + 7) ) ,
			     SSE2_MULCONJ( s1 , *(pw + 4) ) ) ;
  __m128d tmp5 = _mm_sub_pd( SSE2_MULCONJ( s0 , *(pw + 8) ) ,
			     SSE2_MULCONJ( s1 , *(pw + 5) ) ) ;	     
  *( pw + 3 ) = tmp0 ;
  *( pw + 4 ) = tmp1 ;
  *( pw + 5 ) = tmp2 ;
  *( pw + 6 ) = tmp3 ;
  *( pw + 7 ) = tmp4 ;
  *( pw + 8 ) = tmp5 ;
  // compute the new matrix "U -> U * su(2)_1" these are columnal!
  tmp0 = _mm_add_pd( SSE2_MULCONJ( s0 , *(pU + 1) ) ,
		     SSE2_MULCONJ( s1 , *(pU + 2) ) ) ;
  tmp1 = _mm_sub_pd( SSE2_MUL( s0 , *(pU + 2) ) ,
		     SSE2_MUL( s1 , *(pU + 1) ) ) ;
  tmp2 = _mm_add_pd( SSE2_MULCONJ( s0 , *(pU + 4) ) ,
		     SSE2_MULCONJ( s1 , *(pU + 5) ) ) ;
  tmp3 = _mm_sub_pd( SSE2_MUL( s0 , *(pU + 5) ) ,
		     SSE2_MUL( s1 , *(pU + 4) ) ) ;
  tmp4 = _mm_add_pd( SSE2_MULCONJ( s0 , *(pU + 7) ) ,
		     SSE2_MULCONJ( s1 , *(pU + 8) ) ) ;
  tmp5 = _mm_sub_pd( SSE2_MUL( s0 , *(pU + 8) ) ,
		     SSE2_MUL( s1 , *(pU + 7) ) ) ;
  *( pU + 1 ) = tmp0 ;
  *( pU + 2 ) = tmp1 ;
  *( pU + 4 ) = tmp2 ;
  *( pU + 5 ) = tmp3 ;
  *( pU + 7 ) = tmp4 ;
  *( pU + 8 ) = tmp5 ;
  return ;
}

static inline void
rotation3( GLU_complex *U , 
	   GLU_complex *w )
{
  __m128d *pU = (__m128d*)U ;
  __m128d *pw = (__m128d*)w ;
  register __m128d s0 = _mm_add_pd( *(pU + 8) , SSE2_CONJ( *(pU + 0) ) ) ;
  register __m128d s1 = _mm_sub_pd( *(pU + 6) , SSE2_CONJ( *(pU + 2) ) ) ;
  register __m128d sum = _mm_add_pd( _mm_mul_pd( s0 , s0 ) ,
				     _mm_mul_pd( s1 , s1 ) ) ;
  sum = _mm_add_pd( sum , _mm_shuffle_pd( sum , sum , 1 ) ) ;
  sum = _mm_sqrt_pd( sum ) ;
  s0  = _mm_div_pd( s0 , sum ) ;
  s1  = _mm_div_pd( s1 , sum ) ;
  __m128d tmp0 = _mm_sub_pd( SSE2_MULCONJ( s0 , *(pw + 0) ) ,
			     SSE2_MULCONJ( s1 , *(pw + 6) ) ) ;
  __m128d tmp1 = _mm_sub_pd( SSE2_MULCONJ( s0 , *(pw + 1) ) ,
			     SSE2_MULCONJ( s1 , *(pw + 7) ) ) ;
  __m128d tmp2 = _mm_sub_pd( SSE2_MULCONJ( s0 , *(pw + 2) ) ,
			     SSE2_MULCONJ( s1 , *(pw + 8) ) ) ;	
  __m128d tmp3 = _mm_add_pd( SSE2_MUL( s0 , *(pw + 6) ) ,
			     SSE2_MUL( s1 , *(pw + 0) ) ) ;
  __m128d tmp4 = _mm_add_pd( SSE2_MUL( s0 , *(pw + 7) ) ,
			     SSE2_MUL( s1 , *(pw + 1) ) ) ;
  __m128d tmp5 = _mm_add_pd( SSE2_MUL( s0 , *(pw + 8) ) ,
			     SSE2_MUL( s1 , *(pw + 2) ) ) ;
  *( pw + 0 ) = tmp0 ;
  *( pw + 1 ) = tmp1 ;
  *( pw + 2 ) = tmp2 ;
  *( pw + 6 ) = tmp3 ;
  *( pw + 7 ) = tmp4 ;
  *( pw + 8 ) = tmp5 ;
  // compute the new matrix "U -> U * su(2)_1" these are columnal!
  tmp0 = _mm_sub_pd( SSE2_MUL( s0 , *(pU + 0) ) ,
		     SSE2_MUL( s1 , *(pU + 2) ) ) ;
  tmp1 = _mm_add_pd( SSE2_MULCONJ( s0 , *(pU + 2) ) ,
		     SSE2_MULCONJ( s1 , *(pU + 0) ) ) ;
  tmp2 = _mm_sub_pd( SSE2_MUL( s0 , *(pU + 3) ) ,
		     SSE2_MUL( s1 , *(pU + 5) ) ) ;
  tmp3 = _mm_add_pd( SSE2_MULCONJ( s0 , *(pU + 5) ) ,
		     SSE2_MULCONJ( s1 , *(pU + 3) ) ) ;
  tmp4 = _mm_sub_pd( SSE2_MUL( s0 , *(pU + 6) ) ,
		     SSE2_MUL( s1 , *(pU + 8) ) ) ;
  tmp5 = _mm_add_pd( SSE2_MULCONJ( s0 , *(pU + 8) ) ,
		     SSE2_MULCONJ( s1 , *(pU + 6) ) ) ;
  *( pU + 0 ) = tmp0 ;
  *( pU + 2 ) = tmp1 ;
  *( pU + 3 ) = tmp2 ;
  *( pU + 5 ) = tmp3 ;
  *( pU + 6 ) = tmp4 ;
  *( pU + 8 ) = tmp5 ;
  return ;
}

#else

//  The first of the su(2) subgroups for su(3) is
//
//  |  s0   s1  0 |
//  | -s1*  s0* 0 |
//  |  0     0  1 |
//
static void
rotation1( GLU_complex U[ NCNC ] , 
	   GLU_complex w[ NCNC ] )
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
static void
rotation2( GLU_complex U[ NCNC ] ,
	   GLU_complex w[ NCNC ] )
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
static void
rotation3( GLU_complex U[ NCNC ] ,
	   GLU_complex w[ NCNC ] )
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
#endif

#elif NC > 3

// NC generic givens rotations
static void
rotation( GLU_complex U[ NCNC ] , 
	  GLU_complex w[ NCNC ] ,
	  const size_t su2_index )
{
  const size_t a = Latt.su2_data[ su2_index ].idx_a ;
  const size_t b = Latt.su2_data[ su2_index ].idx_b ;
  const size_t c = Latt.su2_data[ su2_index ].idx_c ;
  const size_t d = Latt.su2_data[ su2_index ].idx_d ;

  register GLU_complex s0 = U[a] + conj( U[d] ) ;
  register GLU_complex s1 = U[b] - conj( U[c] ) ;

  // I know how to hypot this, but don't know whether it is worth it
  const double scale = 1.0 / sqrt( creal(s0)*creal(s0) + cimag(s0)*cimag(s0) + \
				   creal(s1)*creal(s1) + cimag(s1)*cimag(s1) ) ;
  s0 *= scale ;
  s1 *= scale ;

  // these overwrite w and U, utilising the fact that most 
  // of the rotation matrix is filled with zeros, ONLY does 
  // the necessary operations.
  // Factor of 10x for SU(8) for doing this, old code is back 
  // in the SVN for checking
  shortened_su2_multiply( w , s0 , s1 , su2_index ) ;
  shortened_su2_multiply_dag( U , s0 , s1 , su2_index ) ;

  return ;
}
#endif

static inline GLU_bool
link_is_zero( const GLU_complex U[ NCNC ] )
{
  register double res = 0.0 ;
#if (NC == 3)
  res += creal(U[0])*creal(U[0]) + cimag(U[0])*cimag(U[0]);
  res += creal(U[1])*creal(U[1]) + cimag(U[1])*cimag(U[1]);
  res += creal(U[2])*creal(U[2]) + cimag(U[2])*cimag(U[2]);
  // middle
  res += creal(U[3])*creal(U[3]) + cimag(U[3])*cimag(U[3]);
  res += creal(U[4])*creal(U[4]) + cimag(U[4])*cimag(U[4]);
  res += creal(U[5])*creal(U[5]) + cimag(U[5])*cimag(U[5]);
  // final
  res += creal(U[6])*creal(U[6]) + cimag(U[6])*cimag(U[6]);
  res += creal(U[7])*creal(U[7]) + cimag(U[7])*cimag(U[7]);
  res += creal(U[8])*creal(U[8]) + cimag(U[8])*cimag(U[8]);
#else
  for( int i = 0 ; i < NCNC ; i++ ) {
    res += creal(U[i])*creal(U[i]) + cimag(U[i])*cimag(U[i]);
  }
#endif
  if( res < PREC_TOL ) {
    return GLU_TRUE ;
  }
  return GLU_FALSE ;
}

// Cabibbo-Marinari projections for trace maximisation
void
givens_reunit( GLU_complex U[ NCNC ] ) 
{
  GLU_complex w[ NCNC ] GLUalign ;
  double trace_new , trace_old ;

  #ifdef OBC_HACK
  if( link_is_zero( U ) == GLU_TRUE ) {
    return ;
  }
  #endif

  #ifdef GIVE_PRECOND
  // MILC precondition,
  // provides a starting guess; the reunitarised path. This is seen 
  // to consistently reduce the number of iterations needed.
  // Seems to work (i.e. reduce breaking of gauge invariance
  // for APE smearing for large NC) Big slow down though due to gram 
  // schmidt and nested SU(2) updates; for SU(3) looks pretty good though 
  // suggest using n-ape for large NC though, is much stabler
  GLU_complex temp[ NCNC ] GLUalign ;
  equiv( temp , U ) ;
  equiv( w , U ) ;
  gram_reunit( w ) ;
  multab_dag( U , temp , w ) ;
  #else
  identity( w ) ;
  #endif

  speed_trace_Re( &trace_old , U ) ;

  size_t i ;
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

    GLU_complex rot[ NCNC ] GLUalign , tmp[ NCNC ] GLUalign ;
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
    // su2 only ever needs one update ...                                       
    break ;

    #else

    // NC-generic cabbibo-marinari update over the su2 subgroups
    size_t su2_index ;
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

// clean up for local scope
#ifdef GIVE_PRECOND
  #undef GIVE_PRECOND
#endif
