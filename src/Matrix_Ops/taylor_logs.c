/*
Copyright 2013-2025 Renwick James Hudspith

    This file (taylor_logs.c) is part of GLU.

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
   @file taylor_logs.c
   @brief logarithms routines and methods to compute the lie-fields of our link matrices

   these are all pretty bad and I would steer clear of them if at all possible
 */
#include "Mainfile.h"

#include "invert.h"

#define ADHOC_PRINCIPAL

/**
   @param NMAX
   @brief maximum number of log iterations
 */
static const size_t NMAX = 25 ; // set this quite high to give it a chance

/**
   @param TAYLOR_SET
   @brief have we set the Taylor expansion coefficients?
 */
static int TAYLOR_SET ;

/**
   @param fact
   @breif Taylor expansion coefficients precomputes
 */
static double *fact ;

// Computes the square-root of the matrix A
//  using a denman-beavers iteration routine
//  A must be invertible for this to make any sense
static int
denman_rootY( GLU_complex Y[ NCNC ] )
{
  GLU_complex M[ NCNC ] , INVM[ NCNC ] ;
  // lets see if the other one makes more sense
  GLU_real tol = 1.0 ;
  size_t i , iters = 0 ;
  for( i = 0 ; i < NCNC ; i++ ) { M[i] = Y[i] ; }

  while( tol > PREC_TOL && iters < 25 ) { // aggressive

    // invert M into INVM
    inverse( INVM , M ) ;

    // M^-1 -> ( 1 + M^-1 ) :: M -> 1 + M
    for( i = 0 ; i < NC ; i++ ) {
      INVM[ i*(NC+1) ] += 1.0 ;
      M[ i*(NC+1) ]    += 1.0 ;
    }
    multab_atomic_right( Y , INVM ) ;

    // Y <- 0.5 Y ( 1 + M^-1 )
    // M <- 0.25 ( M + 2I + M^-1 )
    tol = 0. ;
    for( i = 0 ; i < NCNC ; i++ ) {
      Y[i] *= 0.5 ;
      M[i] = 0.25 * ( M[i] + INVM[i] ) ;
    }
    register GLU_complex c ;
    for( i = 0 ; i < NC ; i++ ) {
      c = M[ i*(NC+1) ] ;
      tol += creal( c ) * creal( c ) + cimag( c ) * cimag( c ) ;
    }
    tol = fabs( tol - NC ) ;

    // should complain here
    if( iters == 24 ) {
      fprintf( stderr , "[TAYLOR_LOGS] Denman root Y ill convergence %e \n" , 
	       tol ) ;
      return GLU_FAILURE ;
    }

    iters++ ;
  }
  // that works most of the time
  return GLU_SUCCESS ;
}

// Computes the inverse square-root of the matrix A
// using a denman-beavers iteration routine.
// Again, A must be invertible. Passes by reference both
// the inverse square root "Z" and the square root "Y"
static int
denman_rootZ( GLU_complex *__restrict Z )
{
  GLU_complex M[ NCNC ] , INVM[ NCNC ] ;
  GLU_real tol = 1.0 ;
  size_t i ,iters = 0 ;
  for( i = 0 ; i < NCNC ; i++ ) {
    M[i] = Z[i] ;
    Z[i] = ( i%(NC+1) ) ? 0.0 : 1.0 ;
  }
  // computes the inverse square root of A (Z), and the square root Y.
  while( tol > PREC_TOL && iters < 25 ) {

    // invert M into INVM
    inverse( INVM , M ) ;

    // add the identity
    for( i = 0 ; i < NC ; i++ ) {
      INVM[ i*(NC+1) ] += 1.0 ;
      M[ i*(NC+1) ]    += 1.0 ;
    }
    multab_atomic_left( Z , INVM ) ;

    // Z <- 0.5 ( 1 + M^-1 ) Z
    // M <- 0.25 ( M + 2I + M^-1 )
    for( i = 0 ; i < NCNC ; i++ ) {
      Z[i] *= 0.5 ; 
      M[i] = 0.25 * ( M[i] + INVM[i] ) ;
    }
    // make sure it is positive definite
    register GLU_complex c ;
    tol = 0. ;
    for( i = 0 ; i < NC ; i++ ) {
      c = M[ i*(NC+1) ] ;
      tol += creal( c ) * creal( c ) + cimag( c ) * cimag( c ) ;
    }
    tol = fabs( tol - NC ) ;

    // should complain here
    if( iters == 24 ) {
      fprintf( stderr , "[TAYLOR_LOGS] Denman root Z ill convergence %e \n" , 
	       tol ) ;
      return GLU_FAILURE ;
    }

    iters++ ;
  }
  // that works!
  return GLU_SUCCESS ;
}

// precomputation, this should be put in Mainfile.c
static void
precompute_taylor_factors( void ) 
{
  // make sure this isn't getting overwritten by different threads
#pragma omp single
  {
    fact = malloc( NMAX * sizeof( double ) ) ;
    size_t n ;
    for( n = 0 ; n < NMAX ; n++ ) {
      fact[ n ] = 1.0 / (double)( 2 * ( n ) - 1 ) ;
    }
    TAYLOR_SET = GLU_TRUE ;
  }
  return ;
}

static int
log_asinh( GLU_complex *__restrict Q , 
	   const GLU_complex *__restrict U ,
	   const size_t NROOTS )
{
  // should I allocate these?, possibly
  GLU_complex y[ NCNC ] , z[ NCNC ] , zz[ NCNC ] ;
  size_t roots , i ;

  // y is our workspace
  equiv( y , U ) ;

  // do 3 some rootings
  for( roots = 0 ; roots < NROOTS ; roots++ ) {
    if( denman_rootY( y ) == GLU_FAILURE ) {
      return GLU_FAILURE ;
    }
  }
  // ok, now compute z = ( U - U^{\dagger} ) / 2
  AntiHermitian_proj( z , y ) ;
  // compute z^2 is the true expansion parameter
  multab( zz , z , z ) ;

  // use pade approximation for asinh
  GLU_complex numerator[ NCNC ] , denominator[ NCNC ] ;

  // asinh using the 6,6 pade
  const int NPADE = 6 ;
  const double num[ 6 ] = { -0.166666666666667 ,  
			    -0.415188764046858 , 
			    -0.371823456681141 ,
			    -0.143528120068249 ,
			    -0.022147200233084 ,
			    -0.000924315649343 } ;
  const double dum[ 6 ] = { 2.941132584281149 ,
			    3.286593260156224 ,
			    1.734623983356869 ,
			    0.435037519640900 ,
			    0.045119586710825 ,
			    0.001245902971704 } ;

  // compute the last term
  for( i = 0 ; i < NCNC ; i++ ) {
    numerator[ i ]   = num[ NPADE-1 ] * zz[ i ] ;
    denominator[ i ] = dum[ NPADE-1 ] * zz[ i ] ;
  }
  // horner's rule for the series of numerator 
  // and denominator
  for( roots = NPADE-1 ; roots == 0 ; roots-- ) {
    for( i = 0 ; i < NC ; i++ ) {
      numerator[i*(NC+1)]   += num[roots] ;
      denominator[i*(NC+1)] += dum[roots] ;
    }
    multab_atomic_left( numerator , zz ) ;
    multab_atomic_left( denominator , zz ) ;
  }
  // add one to the denominator
  for( i = 0 ; i < NC ; i++ ) {
    denominator[i*(NC+1)] += 1.0 ;
  }
  inverse( y , denominator ) ;

  multab_atomic_right( numerator , y ) ;

  // add one to it
  for( i = 0 ; i < NC ; i++ ) {
    numerator[i*(NC+1)] += 1.0 ;
  }

  // multiply on the left by z
  multab_atomic_left( numerator , z ) ;

  // multiply by 2^NROOTS-1
  const GLU_complex SCALE = -I * pow( 2.0 , ( NROOTS ) ) ;
  for( i = 0 ; i < NCNC ; i++ ) {
    Q[i] = SCALE * numerator[i] ;
  }

  // enforce hermiticity to clean up the procedure
  for( i = 0 ; i < NC-1 ; i++ ) { 
    size_t j ;
    for( j = i+1 ; j < NC ; j++ ) {
      Q[ j + i*NC ] = 0.5 * ( Q[ j + i*NC ] + conj( Q[ i + j*NC ] ) ) ;
      Q[ i + j*NC ] = conj( Q[ j + i*NC ] ) ;
    }
  }
  return GLU_SUCCESS ;
}

// logarithm of a unitary matrix by asinh
int
asinh_log( GLU_complex *__restrict Q , 
	   const GLU_complex *__restrict U )
{
  double eps = 1.0 ;
  size_t iters = 0 ;

  GLU_complex Qtemp[ NCNC ] ;
  size_t nroots = 3 ;

  if( log_asinh( Qtemp , U , nroots ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }

  while( eps > 1E-5 && iters < 32 ) {
    nroots+=4 ;

    if( log_asinh( Q , U , nroots ) == GLU_FAILURE ) {
      return GLU_FAILURE ;
    }

    size_t i ;
    eps = 0.0 ;
    for( i = 0 ; i < NCNC ; i++ ) {
      eps += cabs( Qtemp[i] - Q[i] ) ;
    }
    eps /= (NCNC) ;
    iters++ ;
    equiv( Qtemp , Q ) ;
  }
  return GLU_SUCCESS ;
}

//  Analytic brute force Logarithm, it is v. expensive!
//
//  It computes 
//
//  log( U^{1/2^nroots} ) 2^{nroots} = Q
//
//  The matrix square roots are taken by successive Denman-Beavers iterations. It looks
//  like three of these is the sweet spot so I just keep it there.
//
//  Ok, then we do the usually log-chicanery
//
//  log( U ) = y = log( U - 1 ) / log( U + 1 ) series expansion, which is
//
//                            n
//                           ---
//								
//  log( U ) =  y * ( 1.0  +  -  y^(2i-1) * 1/( 2 * i - 1 ) )
//                           /
//                           ---
//                           i=2
//
//
//   Again, the matrix U must be invertible for this to make sense as we use the
//   inverse to accelerate the series. As a note, I do the horner's rule correctly
//   i.e. backwards, with the smallest terms added first to account for roundoff
//   in the matrix power series.
int
brute_force_log( GLU_complex *Q , 
		 const GLU_complex *U ,
		 const size_t NROOTS )
{
  // check to see if we can do precomputations
  if( TAYLOR_SET == GLU_FALSE ) { precompute_taylor_factors( ) ; } 

#if NC < 10
  const double SCALING = pow( 2.0 , ( NROOTS ) ) ;
#else
  const double SCALING = pow( 2.0 , ( 10 ) ) ;
#endif

  // takes nested square roots and puts in the variable "y"
  // should I allocate these?, possibly
  GLU_complex y[ NCNC ] , yy[ NCNC ] , pert[ NCNC ] ;
  equiv( y , U ) ;

  size_t roots , i ;
  // do some rootings ... bit
#if NC < 10
  for( roots = 0 ; roots < NROOTS ; roots++ ) {
#else
  for( roots = 0 ; roots < 10 ; roots++ ) {
#endif
    if( denman_rootY( y ) == GLU_FAILURE ) {
      return GLU_FAILURE ;
    }    

    // durr says that all the eigenvalues should be real, 
    // don't have an NC-generic eigensolver yet
    #ifdef CHECK_ROOTS
    GLU_complex Z[ NCNC ] ;
    Eigenvalues( Z , y ) ;
    int flag = 0 ;
    for( i = 0 ; i < NC ; i++ ) {
      if( creal( Z[i] ) < 0.0 ) { flag = 1 ; break ; } 
    }
    if( flag == 0 ) { break ; } 
    #endif
  }
 
  // yy = U^{1/2^nroots} - 1
  //  y = U^{1/2^nroots} + 1
  for( i = 0 ; i < NCNC ; i++ ) {
    yy[i] = y[i] ;
    if( i%(NC+1) == 0 ) {
      yy[i] -= 1.0 ;
      y[i]  += 1.0 ;
    }
  }
  // invert x + 1
  inverse( Q , y ) ;

  // recall :: y is the n^th square root of U + I
  //        :: yy is the n^th square root of U - I
  //	    :: QP is the inverse of y
  //
  //  y becomes the matrix form of (1-x)/(1+x) where x is the n^th square root
  //  yy is the square of this
  //  we then multiply by 2 to absorb the leading factor in the series expansion
  multab( y , yy , Q ) ;  // compute y
  multab( yy , y , y ) ;  // compute y^2
  for( i = 0 ; i < NCNC ; i++ ) { y[i] *= 2.0 ; }

  double err = 1.0 ;
  size_t n = 4 , j ;
  // use an estimate for the convergence as the last term in the series
  while( n <= NMAX ) {
    // power of the yy matrix
    matrix_power( pert , yy , n + 1 ) ;
    for( j = 0 ; j < NCNC ; j++ ) { 
      pert[ j ] *= fact[ n + 1 ] ;
    }
    // multiply finally by the parameter "y"
    multab_atomic_left( pert , y ) ;
    // compute the error between the two series evaluations
    err = 0.0 ;
    for( j = 0 ; j < NCNC ; j++ ) { 
      err += cabs( pert[j] ) ; 
    }
    err *= SCALING ;
    // this whole code should probably return failure
    if( n >= NMAX ) {
      fprintf( stderr , "[TAYLOR LOG] Not converging :: %e \n" , err ) ;
      return GLU_FAILURE ;
    }
    // increment the series by two to avoid doing too many matrix multiplies. Tune?
    if( err < 0.001*PREC_TOL ) break ;
    n ++ ;
  }

  // compute the next order term in the series in step Q is the n+1 series approx
  //double fact = 1.0 / ( 2. * (n+1) - 1.0 ) ;
  for( j = 0 ; j < NCNC ; j++ ) { 
    Q[j] = yy[j] * fact[ n+1 ] ;
  }

  // Now we have an estimate for when the last term in the series is negligible
  
  // usual horner's rule
  for ( i = n ; i > 1 ; i-- ) {
    //fact = 1.0 / ( 2. * i - 1.0 ) ;
    for( j = 0 ; j < NC ; j++ ) { 
      Q[j*(NC+1)]  += fact[ i ] ; 
    }
    //xx = yy * ( fact + xx ) ;
    multab_atomic_left( Q , yy ) ;
  }
  // xx = y * ( 1.0 + xx ) is the matrix iQ / ( 2^nroots )
  for( j = 0 ; j < NC ; j++ ) { 
    Q[j*(NC+1)] += 1.0 ; 
  }
  
  // multiply finally by the parameter "y"
  multab_atomic_left( Q , y ) ;

  // rescale 2^{NROOTS} outside of the subtraction and ensure tracelessness
  // by enforcing traclessness we break the exponential map slightly
  #ifdef ADHOC_PRINCIPAL
  register const GLU_complex tr = trace ( Q ) / ( NC ) ;
  #endif
  register const GLU_complex f = -I*SCALING ;
  for( i = 0 ; i < NCNC ; i++ ) { 
    #ifdef ADHOC_PRINCIPAL
    if( i%(NC+1)==0 ) { Q[i] -= tr ; }
    #endif
    Q[i] *= f ; 
  }

  // enforce hermiticity to clean up the procedure
  for( i = 0 ; i < NC-1 ; i++ ) { 
    for( j = i+1 ; j < NC ; j++ ) {
      Q[ j + i*NC ] = 0.5 * ( Q[ j + i*NC ] + conj( Q[ i + j*NC ] ) ) ;
      Q[ i + j*NC ] = conj( Q[ j + i*NC ] ) ;
    }
  } 
  return GLU_SUCCESS ;
}

// we free this in Mainfile.c
void
free_taylors( void )
{
  if( TAYLOR_SET == GLU_TRUE ) {
    free( fact ) ;
  }
  return ;
}

// rescaled Denman-Beavers projection, used in the n-APE projection
// which creates a differentiable SU(NC) projection from APE links
int
nape_reunit( GLU_complex *__restrict U )
{
  GLU_complex Z[ NCNC ] ;
  multabdag( Z , U , U ) ;
  // compute determinant shifts
  const GLU_complex dtvvd = cpow( det( Z ) , -1./(double)NC ) ;
  const GLU_complex dtv = cpow( det( U ) , -1./(double)NC ) ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) { 
    Z[i] *= dtvvd ;
    U[i] *= dtv ;
  }
  // Z is the inverse square root of  U.U^{\dagger} 
  if( denman_rootZ( Z ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  // atomically does U -> U Z
  multab_atomic_right( U , Z ) ;

  return GLU_SUCCESS ;
}
