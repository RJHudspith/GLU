/*
Copyright 2013-2025 Renwick James Hudspith

    This file (solver.c) is part of GLU.

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
   @file solver.c
   @brief solves the characteristic equations for SU(2) and SU(3) matrices
 */
#include "Mainfile.h"
#include "invert.h"

#ifdef HAVE_LAPACKE_H
  #include <lapacke.h>
#elif defined HAVE_GSL
  #include <gsl/gsl_complex.h>
  #include <gsl/gsl_complex_math.h>
  #include <gsl/gsl_matrix.h>
  #include <gsl/gsl_eigen.h>
#endif

// I * sqrt(3)
#define rr3 (I * 1.7320508075688772)

#if NC==3
// CALCULATES (one of) THE CUBE ROOT(s) //
static inline void 
cubert( double complex *__restrict res ,
	const double complex z )
{
  register const double real = cabs( z ) ;
  register const double angle = carg( z ) / 3. ; 
  *res = cbrt( real ) * ( cos( angle ) + I * sin( angle ) ) ; 
  return ;
}
#endif

#if !(defined HAVE_GSL_H || defined HAVE_LAPACKE_H) && ( NC > 3 ) 

// deflates a matrix using eigenvector v ... hotelling deflation
static void
deflate( GLU_complex B[ NCNC ] , 
	 const GLU_complex A[ NCNC ] , 
	 const GLU_complex v[ NC ] ,
	 const double evalue )
{
  double norm = 0.0 ;
  size_t j ;
  for( j = 0 ; j < NC ; j++ ) {
    norm += creal( v[ j ] ) * creal( v[ j ] ) +\
      cimag( v[j] ) * cimag( v[j] ) ;
  }
  outerproduct( B , v , v ) ;
  norm = evalue / norm ;
  for( j = 0 ; j < NCNC ; j++ ) {
    B[ j ] = A[ j ] - norm * B[ j ] ;
  }
  return ;
}

/**
   @fn static GLU_complex rayliegh_quotient( GLU_complex v[ NC ] , const GLU_complex A[ NCNC ] )
   @brief computes the eigenvalues of hermitian matrix A
   @param v :: the eigenvectors
   @param A :: the hermitian matrix
   @return the largest eigenvalue of A
   @warning only works for hermitian or symmetric A
 */
static double
rayliegh_quotient( GLU_complex v[ NC ] , 
		   const GLU_complex A[ NCNC ] )
{
  // rayliegh quotient iter ...
  GLU_complex t[ NCNC ] GLUalign , t2[ NCNC ] GLUalign , b[ NC ] ;
  GLU_complex numerator , sum ;
  double evalue = 1.0 , err = 1.0 , norm , new ;
  size_t i , j , iterations = 0 ;
  // close upon convergence
  while( err > PREC_TOL && iterations < 10 ) {
    equiv( t , A ) ;
    // compute t = ( A - evalue*I )
    for( i = 0 ; i < NC ; i++ ) {
      t[ i*(NC+1) ] -= evalue ;
      // set b == v
      b[ i ] = v[ i ] ;
    }
    // invert it into t2
    if( inverse( t2 , t ) == GLU_FAILURE ) break ;
    // multiply by v
    norm = 0.0 ;
    for( i = 0 ; i < NC ; i++ ) {
      for( j = 0 ; j < NC ; j++ ) {
	v[i] += ( t2[ j+i*NC ] * b[j] ) ;
      }
      norm += creal( v[ i ] ) * creal( v[i] ) +	\
	cimag( v[i] ) * cimag( v[i] ) ;
    }
    norm = 1.0 / norm ;
    // set b = v / norm
    for( i = 0 ; i < NC ; i++ ) { v[ i ] *= norm ; }
    // multiply v^{\dagger}Av
    numerator = 0.0 ;
    for( i = 0 ; i < NC ; i++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < NC ; j++ ) {
	sum += A[ j + i*NC ] * v[ j ] ;
      }
      numerator += conj( v[ i ] ) * sum ;
    }
    // compute evalue = v.A.v / v.v ( norm!! )
    new = creal( numerator / norm ) ;
    err = fabs( evalue - new ) ; 
    evalue = new ; // set eigenvalue to the newest computation
    iterations++ ;
  }
  return evalue ;
}
#endif

#if NC==3
// one that returns the root
static inline double complex 
squarert( const double complex z ,
	   const double complex R )
{
  register const double real = cabs( z )  ;
  register const double angle = carg( z ) * 0.5 ;
  const double temp = sqrt( real ) * ( cos( angle ) + I * sin( angle ) ) ; 
  return ( creal( conj( R ) * temp ) < 0 ) ? -temp : temp ;
}
#endif

#if !(defined HAVE_GSL_H || defined HAVE_LAPACKE_H) && ( NC > 3 ) 
/**
   @fn static void rayliegh_evalues( double z[ NC ] , const GLU_complex A[ NCNC ] )
   @brief computes the eigenvalues of A using rayliegh quotient iteration and hotelling deflation
   @param z :: the eigenvalues
   @param A :: some hermitian matrix
 */
static void
rayliegh_evalues( double z[ NC ] ,
		  const GLU_complex A[ NCNC ] )
{
  GLU_complex B[ NCNC ] GLUalign , temp[ NCNC ] GLUalign , v[ NC ] ;
  register double sum = 0.0 , oneONC = 1.0/sqrt((double)NC) ;
  double evalue ;
  size_t i , j ;
  equiv( temp , A ) ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    // initial guess for the eigenvectors
    for( j = 0 ; j < NC ; j++ ) { v[ j ] = oneONC ; }
    // compute the largest eigenvalue by the rayliegh quotient
    evalue = rayliegh_quotient( v , temp ) ; 
    z[i] = (double)evalue ;
    // compute B = A - evalue * ( v v^\dagger) / ( v^\dagger v )
    // is hotelling deflation as I understand it
    deflate( B , temp , v , evalue ) ;
    equiv( temp , B ) ;
    // compute the sum so I don't have to do this for the last one
    sum += z[i] ;
  }
  z[i] = -sum ;
  return ;
}
#endif

// I have taken this method from the appendices of
// arXiv:hep-lat/0409141v2
// computes the evalues of 1-U instead
#if NC == 3
static void
Eigenvalues_su3_stab( double complex z[ NC ] , 
		      const GLU_complex U[ NCNC ] ) 
{ 
  const double complex a = 1.0 - (double complex)trace( U ) / 3. ;
  const double complex b = -conj( a ) ;
  const double complex p = ( b + a * ( 2.0 - a ) ) ;
  const double complex q = ( a * ( 1.0 - b - a * 2.0 * ( 1.0 - a / 3. ) ) + b ) * 1.5 ;
  double complex u1 ;
  double complex res = squarert( q*q + p*p*p , p ) ;
  cubert( &u1 , res - q ) ;
  const double complex u2 = ( creal( u1 ) != 0.0 ) ? -p / ( u1 ) : 0.0 ;
  register const double complex plus = u1 + u2 ;
  register const double complex minus = rr3 * ( u1 - u2 ) ;
  z[0] = 1.0 + plus - a ;
  z[1] = 1.0 + 0.5 * ( minus - plus ) - a ;
  z[2] = z[1] - minus ;
  return ;
}
#endif

// solves the characteristic equation for a general 3*3 matrix. //
void
Eigenvalues( double complex z[ NC ] ,
	     const GLU_complex *__restrict U  )
{
#if NC == 3
  /// gives the roots for x^3+ax^2+bx+c=0 using cardano's method ///
  GLU_complex temp ;
  speed_det( &temp , U ) ; // should exit if the det is 0?
  GLU_complex inv[ NCNC ] GLUalign ; 
  inverse( inv , U ) ; 
  double complex res = (double complex)temp ;
  const double complex trinv = (double complex)trace( inv ) ; 
  const double complex a = -(double complex)trace( U )/3. ; 
  const double complex b = (double complex)res * trinv ; 
  const double complex c = -(double complex)res ; 
  const double complex Q = ( a*a - b/3. ) ; 
  const double complex R = a*a*a - 0.5*a*b + 0.5*c ;
  cubert( &res , R + csqrt( R*R - Q*Q*Q ) ) ; 
  const double complex A = -(double complex)res ; 
  const double complex B = Q / A ; //( A != 0. ? Q/A : 0. ) ; 
  const double complex plus = A + B ;
  const double complex minus = rr3 * ( A - B ) ; 
  z[0] = plus - a ; 
  z[1] = -0.5*plus - a + 0.5*minus ; 
  z[2] = -0.5*plus - a - 0.5*minus ; 
#elif NC == 2
  const double complex a = trace( U ) ;
  const double complex b = det( U ) ;  
  z[0] = 0.5 * ( a + csqrt( a*a - 4.*b ) ) ;
  z[1] = b / z[0] ;
#else
  z[0] = U[0] ;
  // need to think about this one .... hmmm. Might need a library as complex QR
  // appears to be pretty hard to do well
  fprintf( stderr , "[Evalues] sorry not implemented yet .. exiting \n" ) ;
  exit(-1) ;
#endif
  return ;
}

// these are all real can do a lot better than this probably, have a generic
// solver at the bottom
void 
Eigenvalues_hermitian( double z[ NC ] , 
		       const GLU_complex U[ NCNC ] )
{
#if NC == 3
  GLU_complex res ; 
  GLU_real trAA ; 
  trace_prod_herm( &trAA , U ) ; 
  const double b = -0.5 * trAA ; 
  speed_det( &res , U ) ; 
  const double complex c = -(double complex)res ; 
  const double Q = -b/3 ; 
  double complex R = c/2 ; 
  double complex temp = c * c/4. + b * b * b/27. ; 
  temp = R + csqrt( temp ) ; 
  // reuse 'R'
  cubert( &R , temp ) ; 
  register const double complex A = -R ; 
  register const double complex B = Q/A ; 
  register const double complex plus = A + B ;
  register const double complex minus = A - B ; 
  z[0] = (double)creal( plus ) ; 
  z[1] = (double)creal( ( -plus + rr3 * minus ) * 0.5 ) ; 
  z[2] = -z[0] - z[1]  ; 
#elif NC == 2
  z[0] = csqrt( creal( U[0] ) * creal( U[0] ) + creal( U[1] ) * creal( U[1] ) + cimag( U[1] ) * cimag( U[1] ) ) ;
  z[1] = -z[0] ;
#else
  #ifdef HAVE_LAPACKE_H
  // sexy eigenvalues of hermitian matrix thanks to the power of lapack(e)
  double complex a[ NCNC ] ;
  for( i = 0 ; i < NCNC ; i++ ) { a[i] = U[i] ; }
  const int n = NC , lda = NC  ;
  int info = LAPACKE_zheev( LAPACK_ROW_MAJOR , 'N' , 'U' ,
			    n , a , lda , z ) ;
  #elif defined HAVE_GSL
  // grim that we have to do all of these mallocs and frees, gsl's matrix stuff is a mess
  gsl_matrix_complex *mat = gsl_matrix_complex_alloc(NC,NC);
  gsl_vector *eval = gsl_vector_alloc(NC) ;
  gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc(NC) ;
  size_t i , j ;
  register double re , im ;
  // pack the matrices
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      re = (double)creal( U[j+NC*i] ) ;
      im = (double)cimag( U[j+NC*i] ) ; 
      gsl_matrix_complex_set( mat , i , j , gsl_complex_rect( re , im ) ) ;
    }
  } 
  gsl_eigen_herm ( mat , eval, w) ; // calculate the eigenvalues
  for( i = 0 ; i < NC ; i++ ) {
    z[i] = gsl_vector_get (eval, i);
  }
  gsl_vector_free( eval ) ;
  gsl_matrix_complex_free( mat ) ;
  gsl_eigen_herm_free( w ) ;
  #else
  // rayliegh quotient iteration, checked against GSL's evaluation
  // uses hotelling deflation to get the rest, and ensures the eigenvalues
  // sum to 0.0 !
  rayliegh_evalues( z , U ) ;
  #endif
#endif
  return ;
}

// SU(NC) variant eigenvalues
void
Eigenvalues_suNC( double complex z[ NC ] , 
		 const GLU_complex U[ NCNC ] ) 
{
#if NC == 3
  if( 1.0 - ( creal( U[0] ) + creal( U[4] ) + creal( U[8] ) )/3.  < 1E-12 ) { 
    Eigenvalues_su3_stab( z , U ) ; 
    return ;
  }
  const double complex a = -(double complex)( U[0] + U[4] + U[8] ) / 3. ;
  register const double complex conjA = conj( a ) ;
  const double complex Q = ( a * a ) + conjA ;
  const double complex R = a * ( Q + conjA * 0.5 ) - 0.5 ; 
  const double complex temp = R * R - Q * Q * Q ; 
  double complex res = squarert( temp , R ) ;
  cubert( &res , R + res ) ; 
  register const double complex A = -res ; 
  register const double complex B = ( A != 0.0 ? Q/A : 0.0 ) ; //Q/A
  // I don't think we will ever get a div 0 error here
  register const double complex plus = A + B ;
  register const double complex minus = rr3 * ( A - B ) ; 
  z[0] = plus - a ; 
  z[1] = ( minus - plus ) * 0.5 - a ; 
  z[2] = z[1] - minus ; 
#elif NC == 2
  const double reU0 = (double)creal( U[0] ) ;
  z[0] = reU0 + I * sqrt ( 1.0 - reU0 * reU0 ) ;
  z[1] = conj( z[0] ) ;
#else
  // not supported, yet. Wrap to eigenvalues ....
  Eigenvalues( z , U ) ;
#endif
  return ;
}

#undef rr3
