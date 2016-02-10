/*
    Copyright 2013 Renwick James Hudspith

    This file (expMat.c) is part of GLU.

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
   @file expMat.c
   @brief matrix exponentiation routines
 */
#include "Mainfile.h"

#include "effs.h"     // compute the f-constants of Cayley-Hamilton
#include "solver.h"   // characteristic equation eigensolver

// sqrt of three appears everywhere
#define r3 (1.7320508075688772935)

// if we don't have easier options we need some routines defined elsewhere
#if !( defined HAVE_LAPACKE_H || defined HAVE_GSL ) && ( NC > 3 )

// for large NC the taylor expansion exponential is used which uses inverse()
#include "invert.h"
#include "gramschmidt.h"

// use the pade approximation! It is MUCH faster
#define USE_PADE

#ifdef SINGLE_PREC
  #define MAX_FACTORIAL 10
#else
  #define MAX_FACTORIAL 20
#endif

// LUT for the factorials
#ifdef USE_PADE
#define PADE_ORDER (4)
static GLU_real *pade ;
#else
static GLU_real *factorial ;
#endif

// initialise the factorial in Mainfile.c
void
init_factorial( void ) 
{
#ifdef USE_PADE
  // mix of 4,4 pade and 3 rootings is nice combination
  pade = ( GLU_real* ) malloc( (PADE_ORDER+1) * sizeof( GLU_real ) ) ;
  pade[ 0 ] = 1.0 ;
  pade[ 1 ] = 1.0/2.0 ;
  pade[ 2 ] = 3.0/28.0 ;
  pade[ 3 ] = 1.0/84.0 ;
  pade[ 4 ] = 1.0/1680.0 ;
#else
  factorial = ( GLU_real* ) malloc ( MAX_FACTORIAL * sizeof( GLU_real ) ) ;
  factorial[ 0 ] = 1.0 ;
  factorial[ 1 ] = 1.0 ;
  factorial[ 2 ] = 0.5 ;
  factorial[ 3 ] = 1.0/6 ;
  factorial[ 4 ] = 1.0/24 ;
  factorial[ 5 ] = 1.0/120 ;
  factorial[ 6 ] = 1.0/720 ;
  factorial[ 7 ] = 1.0/5040 ;
  factorial[ 8 ] = 1.0/40320 ;
  factorial[ 9 ] = 1.0/362880 ;
  #ifndef SINGLE_PREC
  factorial[ 10 ] = 1.0/39916800 ;
  factorial[ 11 ] = 1.0/479001600 ;
  factorial[ 12 ] = 1.0/6227020800 ;
  factorial[ 13 ] = 1.0/87178291200 ;
  factorial[ 14 ] = 1.0/1307674368000 ;
  factorial[ 15 ] = 1.0/20922789888000 ;
  factorial[ 16 ] = 1.0/355687428096000 ;
  factorial[ 17 ] = 1.0/6402373705728000 ;
  factorial[ 18 ] = 1.0/121645100408832000 ;
  factorial[ 19 ] = 1.0/2432902008176640000 ;
  #endif
#endif
  return ;
}

// free the factorial too 
void
free_factorial( void )
{
#ifdef USE_PADE
  free( pade ) ;
#else
  free( factorial ) ;
#endif
  return ;
}

#ifdef USE_PADE
// n,n pade approximation for the exponential
static void
horners_pade( GLU_complex a[ NCNC ] ,
	      const GLU_complex b[ NCNC ] )
{
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) { a[i] = b[i] * pade[PADE_ORDER] ; } 
  for( i = PADE_ORDER-1 ; i > 0 ; i-- ) {    
    add_constant( a , pade[i] ) ; 
    multab_atomic_left( a , b ) ;    
  }
  for( i = 0 ; i < NC ; i++ ) { a[ i*(NC+1) ] += 1.0 ; }
  return ;
}
#else
// horner's expansion for the exponential
static void
horners_exp( GLU_complex a[ NCNC ] ,
	     const GLU_complex b[ NCNC ] ,
	     const size_t n )
{
  size_t i ;
  // I use max iterations 13 for double prec and 8 for single it was before
  // this selection was borne out of heavy testing I would urge you to keep it
  #ifndef SINGLE_PREC
  for( i = 0 ; i < NCNC ; i++ ) { a[i] = b[i] * factorial[13] ; } 
  for( i = 12 ; i > 0 ; i-- ) {
    add_constant( a , factorial[i] ) ; 
    multab_atomic_left( a , b ) ;    
  }
  #else
  for( i = 0 ; i < NCNC ; i++ ) { a[i] = b[i] * factorial[9] ; } 
  for( i = 8 ; i > 0 ; i-- ) {    
    add_constant( a , factorial[i] ) ; 
    multab_atomic_left( a , b ) ;    
  }
  #endif
  for( i = 0 ; i < NC ; i++ ) { a[ i*(NC+1) ] += 1.0 ; }
  return ;
}
#endif

#endif

// exponentiate exactly a hermitian matrix "Q" into SU(NC) matrix "U"
void
exponentiate( GLU_complex U[ NCNC ] , 
	      const GLU_complex Q[ NCNC ] )
{
#if NC == 3
  GLU_real *qq = ( GLU_real* )Q ;
  const double REQ0 = *( qq + 0 ) ;
  const double REQ1 = *( qq + 2 ) ;
  const double IMQ1 = *( qq + 3 ) ;
  const double REQ2 = *( qq + 4 ) ;
  const double IMQ2 = *( qq + 5 ) ;
  const double REQ4 = *( qq + 8 ) ;
  const double REQ5 = *( qq + 10 ) ;
  const double IMQ5 = *( qq + 11 ) ;
  const double REQ8 = *( qq + 16 ) ;

  // speed this up too (use determinant relation)
  const double c1 = ( REQ0 * REQ0 + REQ0 * REQ4 + REQ4 * REQ4		\
		      + REQ1 * REQ1 + IMQ1 * IMQ1			\
		      + REQ2 * REQ2 + IMQ2 * IMQ2			\
		      + REQ5 * REQ5 + IMQ5 * IMQ5 ) * OneO3 ; 
 
  //Iff c0_max < ( smallest representable double) the matrix Q is zero and its
  //exponential is the identity matrix ..
  if( unlikely( c1 < DBL_MIN ) ) {
    *( U + 0 ) = 1. ; 
    *( U + 1 ) = 0. ; 
    *( U + 2 ) = 0. ; 
    *( U + 3 ) = 0. ; 
    *( U + 4 ) = 1. ; 
    *( U + 5 ) = 0. ; 
    *( U + 6 ) = 0. ; 
    *( U + 7 ) = 0. ; 
    *( U + 8 ) = 1. ; 
    return ;
  }

  // will write this out as it can be done cheaper
  // 1/3 * tr AAA is just det( A ) 
  // Below is a quickened determinant
  double c0 =  REQ0  * ( REQ4 * REQ8				\
			 - REQ5  * REQ5  - IMQ5  * IMQ5  ) ;
  // from the middle
  c0 -= REQ1  * ( REQ1 * REQ8 		\
		  - REQ5  * REQ2  - IMQ5  * IMQ2  ) ;
  c0 += IMQ1  * ( - IMQ1 * REQ8		\
		  + REQ5  * IMQ2  - IMQ5  * REQ2 ) ;
  // final column
  c0 += REQ2  * ( - REQ4  * REQ2			\
		  + REQ1 * REQ5  - IMQ1 * IMQ5  ) ;
  c0 -= IMQ2  * ( REQ4  * IMQ2				\
		  - REQ1 * IMQ5  - IMQ1 * REQ5 ) ;

  // so if c0 is negative we flip the sign ...
  const double flag = c0 < 0 ? -1.0 : 1.0 ;
  c0 *= flag ;
 
  // compute the constants c0_max and the root of c1 ...
  const double rc1 = sqrt( c1 ) ;
  const double c0_max = 2. * rc1 * c1 ; 
  const double theta = acos( c0 / c0_max ) * OneO3 ; 
  const double ctheta = cos( theta ) ;
  register const double u = rc1 * ctheta ; 
  register const double w = r3 * rc1 * sin( theta ) ;
  const double uu = u * u  ,  ww = w * w  ,  cw = cos( w ) ; 
  const double denom = 1.0 / ( 9. * uu - ww ) ;
  const double cu = cos( u ) ;
  const double su = sin( u ) ;
  // and I thought double angle formulas were useless!
  //double complex one , two ;
  const double complex one = cu - I * su ;
  double complex two = conj( one ) ; //cu + I * su ;
  two *= two ;
 
  // taylor expand if getting toward the numerically unstable end
  const double E0 = fabs( w ) < SINTOL ? ( 1 - ww / 6. * ( 1 - ww / 20. * ( 1 - ww / 42. ) ) ) : sin( w ) / w ; 

  double complex f0 = ( uu - ww ) * two + one * ( 8. * uu * cw + 2. * I * u * ( 3. * uu + ww ) * E0 ) ; 
  double complex f1 = 2. * u * two - one * ( 2. * u * cw - I * ( 3. * uu - ww ) * E0 ) ; 
  double complex f2 = two - one * ( cw + 3. * I * u * E0 ) ; 

  f0 = denom * ( creal( f0 ) + I * cimag( f0 ) * flag ) ;
  f1 = denom * ( flag * creal( f1 ) + I * cimag( f1 ) ) ;
  f2 = denom * ( creal( f2 ) + I * cimag( f2 ) * flag ) ;

  // QQ[0].
  const double temp0 = REQ0 * REQ0 + REQ1 * REQ1 +	\
    IMQ1 * IMQ1 + REQ2 * REQ2 + IMQ2 * IMQ2 ;
  // QQ[1]
  const double complex temp1 = -REQ1 * ( REQ8 ) + REQ2 * REQ5 + IMQ2 * IMQ5 
    + I * ( REQ5 * IMQ2 - REQ2 * IMQ5 - IMQ1 * REQ8 ) ;
  // QQ[2]
  const double complex temp2 = REQ1 * REQ5 - IMQ1 * IMQ5 - REQ2 * REQ4 +
    I * ( IMQ1 * REQ5 + IMQ5 * REQ1 - IMQ2 * REQ4 ) ;
  // QQ[4]
  const double temp3 = REQ4 * REQ4 + REQ1 * REQ1	\
    + IMQ1 * IMQ1 + REQ5 * REQ5 + IMQ5 * IMQ5 ;
  // QQ[5]
  const double complex temp4 = REQ1 * REQ2 + IMQ2 * IMQ1 - REQ0 * REQ5 +
    I * ( REQ1 * IMQ2 - REQ2 * IMQ1 - REQ0 * IMQ5 ) ;
  // QQ[8]
  const double temp5 = REQ8 * REQ8 + REQ2 * REQ2 +	\
    IMQ2 * IMQ2 + REQ5 * REQ5 + IMQ5 * IMQ5 ;
  // U = f0I + f1 Q + f2 QQ 
  *( U + 0 ) = f0 + f1 * REQ0 + f2 * temp0  ; 
  *( U + 1 ) = f1 * Q[1] + f2 * temp1 ; 
  *( U + 2 ) = f1 * Q[2] + f2 * temp2 ; 
  //
  *( U + 3 ) = f1 * Q[3] + f2 * conj( temp1 ) ; 
  *( U + 4 ) = f0 + f1 * REQ4 + f2 * temp3 ; 
  *( U + 5 ) = f1 * Q[5] + f2 * temp4 ; 
  //
  *( U + 6 ) = f1 * Q[6] + f2 * conj( temp2 ) ; 
  *( U + 7 ) = f1 * Q[7] + f2 * conj( temp4 ) ; 
  *( U + 8 ) = f0 + f1 * REQ8 + f2 * temp5 ; 
#elif NC == 2
  double f0 , f1 ; // f1 is purely imaginary
  // eigenvalues are pretty simple +/- sqrt( |a|^2 + |b|^2 ) Only need one
  const double z = sqrt( creal( Q[0] ) * creal( Q[0] ) +	\
			 creal( Q[1] ) * creal( Q[1] ) +	\
			 cimag( Q[1] ) * cimag( Q[1] ) ) ;

  // have eigenvalues, now for the "fun" bit.
  f0 = cos( z ) ;
  // taylor expand 
  f1 = fabs ( z ) < SINTOLSU2 ? ( 1 - z / 6. * ( 1 - z / 20. * ( 1 - z / 42. ) ) ) : sin( z ) / z ;

  const double complex f1Q0 = I * f1 * creal( Q[0] ) ;
  const double complex f1Q1 = I * f1 * Q[1] ;
  *( U + 0 ) = (GLU_complex)( f0 + f1Q0 ) ;
  *( U + 1 ) = (GLU_complex)( f1Q1 ) ;
  *( U + 2 ) = (GLU_complex)( -conj( f1Q1 ) ) ;
  *( U + 3 ) = (GLU_complex)( f0 - f1Q0 ) ; 
#else
  // hmmm could be a toughy
  #if ( defined HAVE_LAPACKE_H || defined HAVE_GSL )
  double complex f[ NC ] ;
  double z[ NC ] ;
  Eigenvalues_hermitian( z , Q ) ;
  calculate_effs_VDM_herm( f , z ) ;
  // matrix expansion reversing horner's rule
  int i , j ;
  diag( U , f[ NC - 1 ] ) ; 
  for( i = NC-1 ; i > 0 ; i-- ) {
    multab_atomic_left( U , Q ) ; // left multiply U with Q
    for( j = 0 ; j < NC ; j++ ) { U[ j*(NC+1) ] += f[ i-1 ] ; }
  }
  #else
  // exponentiate routine from stephan durr's paper
  // Performs the nesting
  // U = ( exp{ A / DIV^n ) ) ^ ( DIV * n )
  GLU_complex EOLD[ NCNC ] , SN[ NCNC ] ;
  GLU_complex RN_MIN[ NCNC ] , RN[ NCNC ] ;

  // set to zero
  zero_mat( EOLD ) ;

  // set up the divisor and the minimum 
  double sum ;
  const int DIV = 2 , nmin = 3 ;
  size_t j , n ;

  // use precomputed factorials  
  for( n = nmin ; n < 10 ; n++ ) {

    // compute the multiplicative factor ...
    const GLU_real fact = 1.0 / pow( DIV , n ) ;
    const int iter = pow( DIV , n ) ;

    // and the rational approximations
#ifdef USE_PADE
    for( j = 0 ; j < NCNC ; j++ ) { SN[ j ] = ( I * Q[j] * fact ) ; }
    horners_pade( RN , SN ) ;

    for( j = 0 ; j < NCNC ; j++ ) { SN[ j ] *= -1.0 ; }
    horners_pade( RN_MIN , SN ) ; 
#else
    for( j = 0 ; j < NCNC ; j++ ) { SN[ j ] = ( I * Q[j] * fact ) / 2.0 ; }
    horners_exp( RN , SN , 14 ) ;

    for( j = 0 ; j < NCNC ; j++ ) { SN[ j ] *= -1.0 ; }
    horners_exp( RN_MIN , SN , 14 ) ; 
#endif

    inverse( SN , RN_MIN ) ; // uses our numerical inverse
    multab_atomic_right( RN , SN ) ; // gets the correct rational approx

    // and remove the nested scalings ...
    matrix_power( U , RN , iter ) ; // uses a fast-power like routine 
    
    // for the convergence criteria, I use the absolute difference between
    // evaluations
    sum = 0.0 ;
    for( j = 0 ; j < NCNC ; j++ ) {
      sum += (double)cabs( EOLD[j] - U[j] ) ;
      EOLD[ j ] = U[ j ] ;
    }
    sum /= NCNC ;

    // convergence ....
    if( sum < PREC_TOL ) { break ; } 
    // warning for non-convergence ..
    if( n >= ( MAX_FACTORIAL - 1 ) ) { 
      printf( "[EXPONENTIAL] not converging .. %zu %e \n" , n , sum ) ; 
      break ;
    }
  }
  // gramschmidt orthogonalisation just to make sure
  // has been seen to help preserve gauge invariance of log smearing
  reunit2( U ) ;
  #endif
#endif
  return ;
}

// exactly the same as above ,  just calculates the f-functions in-step instead of requiring them
// takes a shortened, Hermitian Q and gives back the SU(NC) matrix U
void
exponentiate_short( GLU_complex U[ NCNC ] , 
		    const GLU_complex Q[ HERMSIZE ] )
{
#if NC == 3
  GLU_real *qq = ( GLU_real* )Q ;
  const double REQ0 = *( qq + 0 ) ;
  const double REQ1 = *( qq + 2 ) ;
  const double IMQ1 = *( qq + 3 ) ;
  const double REQ2 = *( qq + 4 ) ;
  const double IMQ2 = *( qq + 5 ) ;
  const double REQ4 = *( qq + 6 ) ;
  const double REQ5 = *( qq + 8 ) ;
  const double IMQ5 = *( qq + 9 ) ;
  const double REQ8 = -( REQ0 + REQ4 ) ;
  const double c1 = ( REQ0 * -REQ8 + REQ4 * REQ4			\
		      + REQ1 * REQ1 + IMQ1 * IMQ1			\
		      + REQ2 * REQ2 + IMQ2 * IMQ2			\
		      + REQ5 * REQ5 + IMQ5 * IMQ5 ) * OneO3 ;
  
  //  Iff c0_max < ( smallest representable double) the matrix Q is zero and its
  //  exponential is the identity matrix .
  if( unlikely( c1 < DBL_MIN ) ) {
    *( U + 0 ) = 1. ; 
    *( U + 1 ) = 0. ; 
    *( U + 2 ) = 0. ; 
    //
    *( U + 3 ) = 0. ; 
    *( U + 4 ) = 1. ; 
    *( U + 5 ) = 0. ; 
    //
    *( U + 6 ) = 0. ; 
    *( U + 7 ) = 0. ; 
    *( U + 8 ) = 1. ; 
    return ;
  }

  // 1/3 * tr AAA is just det( A )
  // Below is a quickened determinant
  double c0 =  REQ0  * ( REQ4 * REQ8				\
			 - REQ5  * REQ5  - IMQ5  * IMQ5  ) ;
  // from the middle
  c0 -= REQ1  * ( REQ1 * REQ8 		\
		  - REQ5  * REQ2  - IMQ5  * IMQ2  ) ;
  c0 += IMQ1  * ( - IMQ1 * REQ8		\
		  + REQ5  * IMQ2  - IMQ5  * REQ2 ) ;
  // final column
  c0 += REQ2  * ( - REQ4  * REQ2			\
		  + REQ1 * REQ5  - IMQ1 * IMQ5  ) ;
  c0 -= IMQ2  * ( REQ4  * IMQ2				\
		  - REQ1 * IMQ5  - IMQ1 * REQ5 ) ;

  // so if c0 is negative we flip the sign ...
  const double flag = c0 < 0 ? -1.0 : 1.0 ;
  c0 *= flag ;

  // compute the constants c0_max and the root of c1 ...
  const double rc1 = sqrt( c1 ) ;
  const double c0_max = 2. * rc1 * c1 ; 
  const double theta = acos( c0 / c0_max ) * OneO3 ; 
  const double ctheta = cos( theta ) ;
  register const double u = rc1 * ctheta ; 
  register const double w = r3 * rc1 * sin( theta ) ; 
  const double uu = u * u  ,  ww = w * w  ,  cw = cos( w ) ; 
  const double denom = 1.0 / ( 9. * uu - ww ) ;
  const double cu = cos( u ) ;
  const double su = sin( u ) ;
  const double complex one = cu - I * su ;
  double complex two = conj( one ) ;
  two *= two ;

  // taylor expand if getting toward the numerically unstable end
  const double E0 = fabs( w ) < SINTOL ? ( 1 - ww / 6. * ( 1 - ww / 20. * ( 1 - ww / 42. ) ) ) : sin( w ) / w ; 

  double complex f0 = ( uu - ww ) * two + one * ( 8 * uu * cw + 2 * I * u * ( 3 * uu + ww ) * E0 ) ; 
  double complex f1 = 2. * u * two - one * ( 2. * u * cw - I * ( 3 * uu - ww ) * E0 ) ; 
  double complex f2 = two - one * ( cw + 3 * I * u * E0 ) ; 

  f0 = ( creal( f0 ) + I * cimag( f0 ) * flag ) ;
  f1 = ( flag * creal( f1 ) + I * cimag( f1 ) ) ;
  f2 = ( creal( f2 ) + I * cimag( f2 ) * flag ) ;

  f0 *= denom ;
  f1 *= denom ;
  f2 *= denom ;

  // QQ[0].
  const double temp0 = REQ0 * REQ0 + REQ1 * REQ1 +	\
    IMQ1 * IMQ1 + REQ2 * REQ2 + IMQ2 * IMQ2 ;
  // QQ[1]
  const double complex temp1 = -REQ1 * ( REQ8 ) + REQ2 * REQ5 + IMQ2 * IMQ5 
    + I * ( REQ5 * IMQ2 - REQ2 * IMQ5 - IMQ1 * REQ8 ) ;
  // QQ[2]
  const double complex temp2 = REQ1 * REQ5 - IMQ1 * IMQ5 - REQ2 * REQ4 +
    I * ( IMQ1 * REQ5 + IMQ5 * REQ1 - IMQ2 * REQ4 ) ;
  // QQ[4]
  const double temp3 = REQ4 * REQ4 + REQ1 * REQ1	\
    + IMQ1 * IMQ1 + REQ5 * REQ5 + IMQ5 * IMQ5 ;
  // QQ[5]
  const double complex temp4 = REQ1 * REQ2 + IMQ2 * IMQ1 - REQ0 * REQ5 +
    I * ( REQ1 * IMQ2 - REQ2 * IMQ1 - REQ0 * IMQ5 ) ;
  // QQ[8]
  const double temp5 = REQ8 * REQ8 + REQ2 * REQ2 +	\
    IMQ2 * IMQ2 + REQ5 * REQ5 + IMQ5 * IMQ5 ;

  //can really speed this up
  *( U + 0 ) = f0 + f1 * REQ0 + f2 * temp0  ; 
  *( U + 1 ) = f1 * Q[1] + f2 * temp1 ; 
  *( U + 2 ) = f1 * Q[2] + f2 * temp2 ; 
  //
  *( U + 3 ) = f1 * conj( Q[1] ) + f2 * conj( temp1 ) ; 
  *( U + 4 ) = f0 + f1 * REQ4 + f2 * temp3 ; 
  *( U + 5 ) = f1 * Q[4] + f2 * temp4 ; 
  //
  *( U + 6 ) = f1 * conj( Q[2] ) + f2 * conj( temp2 ) ; 
  *( U + 7 ) = f1 * conj( Q[4] ) + f2 * conj( temp4 ) ; 
  *( U + 8 ) = f0 + f1 * REQ8 + f2 * temp5 ; 
#elif NC == 2
  double f0 , f1 ; // f1 is purely imaginary
  // eigenvalues are pretty simple +/- sqrt( |a|^2 + |b|^2 ) Only need one
  const double z = sqrt( creal( Q[0] ) * creal( Q[0] ) +	\
			 creal( Q[1] ) * creal( Q[1] ) +	\
			 cimag( Q[1] ) * cimag( Q[1] ) ) ;

  // have eigenvalues, now for the "fun" bit.
  f0 = cos( z ) ;
  // taylor expand 
  f1 = fabs ( z ) < SINTOLSU2 ? ( 1.0 - z / 6. * ( 1 - z / 20. * ( 1 - z / 42. ) ) ) : sin( z ) / z ;

  const double complex f1Q0 = I * f1 * Q[0] ;
  const double complex f1Q1 = I * f1 * Q[1] ;
  *( U + 0 ) = (GLU_complex)( f0 + f1Q0 ) ;
  *( U + 1 ) = (GLU_complex)( f1Q1 ) ;
  *( U + 2 ) = (GLU_complex)( -conj( f1Q1 ) ) ;
  *( U + 3 ) = (GLU_complex)( f0 - f1Q0 ) ; 
#else
  // hmmm could be a toughy, wrap to exponentiate
  GLU_complex temp[ NCNC ] ;
  rebuild_hermitian( temp , Q ) ;
  exponentiate( U , temp ) ;
#endif
  return ;
}

// clean up r3
#undef r3

// and clear it up
#if !( defined HAVE_LAPACKE_H || defined HAVE_GSL ) && ( NC > 3 )
 #ifdef USE_PADE
  #undef USE_PADE
 #endif
 #undef MAX_FACTORIAL
#endif
