/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (GLU_splines.c) is part of GLU.

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
   @file GLU_splines.c
   @brief cubic spline calculator
 */
#include "Mainfile.h"

#define sgn( a ) ( a < 0.0 ? -1 : 1 ) 

// finds the upper index of target
static size_t
find_idx( const double target , 
	  const double *__restrict X , 
	  const size_t high , 
	  const size_t low )
{
  return ( ( high - low ) < 2 ) ? high :				\
    X[ ( high + low ) >> 1 ] < target ?					\
    find_idx( target , X , high , ( high + low ) >> 1 ) :		\
    find_idx( target , X , ( high + low ) >> 1 , low ) ;
}

// lowest order derivative I support
static void
h2_derivative( double *__restrict der ,
	       const double *__restrict x ,
	       const double *__restrict y ,
	       const size_t N )
{
  // forwards iteration
  // expansion of f(x+a) + f(x+b)
  const double a0 = x[2] - x[0] ;
  const double b0 = x[1] - x[0] ;
  const double num0 = b0*b0*y[2] - a0*a0*y[1] + y[0]*(a0*a0-b0*b0) ;
  der[0] = num0 / ( a0*b0*( b0-a0 ) ) ;

  // and the middle ones
  size_t i ;
  for( i = 1 ; i < N-1 ; i++ ) {
    const double a = x[i+1] - x[i] ;
    const double b = x[i] - x[i-1] ;
    const double num = b*b*y[i+1] - a*a*y[i-1] + y[i]*(a*a-b*b) ;
    der[i] = num / ( a*b*( b+a ) ) ;
  }
  // backwards iteration
  // expansion of f(x-a) + f(x-b)
  const double a1 = x[i] - x[i-2] ;
  const double b1 = x[i] - x[i-1] ;
  const double num1 = b1*b1*y[i-2] - a1*a1*y[i-1] + y[i]*(a1*a1-b1*b1) ;
  der[i] = num1 / ( a1*b1*( a1-b1 ) ) ;
  return ;
}

// lower order derivative
static void
h3_derivative( double *__restrict der ,
	       const double *__restrict x ,
	       const double *__restrict y ,
	       const size_t N )
{
  {
    // expansion of f(x+a) + f(x+b) + f(x+c)
    const double a = x[1] - x[0] ;
    const double b = x[2] - x[0] ;
    const double c = x[3] - x[0] ;
    register double ab = a*b ;
    register double ac = a*c ;
    register double bc = b*c ;
    const double num = bc*bc*(c-b)*y[1] - ac*ac*(-a+c)*y[2] + ab*ab*(-a+b)*y[3] ;
    der[0] = num / ( -ab*(-ac+bc)*(b-c)*(-a+c) ) + (-1./a - 1./b - 1./c ) * y[0] ; 
  }
  // expansion of f(x-a) + f(x+b) + f(x+c)
  size_t i ;
  for( i = 1 ; i < N-2 ; i++ ) {
    const double a = x[i] - x[i-1] ;
    const double b = x[i+1] - x[i] ;
    const double c = x[i+2] - x[i] ;
    register double ab = a*b ;
    register double ac = a*c ;
    register double bc = b*c ;
    const double num = bc*bc*(c-b)*y[i-1] - ac*ac*(a+c)*y[i+1] + ab*ab*(a+b)*y[i+2] ;
    der[i] = num / ( a*b*c*(a+b)*(b-c)*(a+c) ) + (1./a - 1./b - 1./c ) * y[i] ;
  }
  // expansion of f(x-a) + f(x-b) + f(x+c)
  {
    const double a = x[N-2] - x[N-4] ;
    const double b = x[N-2] - x[N-3] ;
    const double c = x[N-1] - x[N-2] ;
    register double ab = a*b ;
    register double ac = a*c ;
    register double bc = b*c ;
    double num = bc*bc*(b+c)*y[N-4] - ac*ac*(a+c)*y[N-3] + ab*ab*(a-b)*y[N-1] ;
    der[N-2] = num / ( ab*c*(a-b)*(b+c)*(a+c) ) + (1./a + 1./b - 1./c)*y[N-2] ;
  }
  // expansion of f(x-a) + f(x-b) + f(x-c)
  {
    const double a = x[N-1] - x[N-4] ;
    const double b = x[N-1] - x[N-3] ;
    const double c = x[N-1] - x[N-2] ;
    register double ab = a*b ;
    register double ac = a*c ;
    register double bc = b*c ;
    const double num = bc*bc*(b-c)*y[N-4] - ac*ac*(a-c)*y[N-3] + ab*ab*(a-b)*y[N-2] ;
    der[N-1] = num / ( -ab*c*(a-b)*(b-c)*(a-c) ) + (1./a + 1./b + 1./c)*y[N-1] ;
  }
  return ;
}

// compute the derivative for the spline
void
spline_derivative( double *__restrict der ,
		   const double *__restrict x ,
		   const double *__restrict y ,
		   const size_t N )
{
  // for these derivative functions to work, need at least 5 points
  if( N < 5 ) { 
    if( N == 3 ) { 
      h2_derivative( der , x , y , N ) ; 
      return ;
    } else if( N == 4 ) { 
      h3_derivative( der , x , y , N ) ; 
      return ;
    }
    printf( "[SPLINE] Too few points to compute a good cubic spline \n" ) ;
    printf( "[SPLINE] Not doing anything \n" ) ;
    return ;
  }

  // +2 forward-looking O(h^4)
  {
    const double a = x[1] - x[0] ;
    const double b = x[2] - x[0] ;
    const double c = x[3] - x[0] ;
    const double d = x[4] - x[0] ;
    const double A = (b*c*d)*(b*c*d)*(-b+c)*(c-d)*(-b+d)*y[1] ;
    const double B = -(a*c*d)*(a*c*d)*(-a+c)*(c-d)*(-a+d)*y[2] ;
    const double C = -(a*b*d)*(a*b*d)*(-a+b)*(-a+d)*(-b+d)*y[3] ;
    const double D = (a*b*c)*(a*b*c)*(-a+b)*(-a+c)*(-b+c)*y[4] ;
    der[0] = ( A+B+C+D )/(-a*-b*c*d*(-a+b)*(-a+c)*(-b+c)*(c-d)*(-a+d)*(-b+d)) +	\
      ( 1./-a + 1./-b - 1./c - 1./d ) * y[0] ;
  }
  // +1 forward-looking O(h^4)
  {
    const double a = x[1] - x[0] ;
    const double b = x[2] - x[1] ;
    const double c = x[3] - x[1] ;
    const double d = x[4] - x[1] ;
    const double A = (b*c*d)*(b*c*d)*(-b+c)*(c-d)*(-b+d)*y[0] ;
    const double B = -(a*c*d)*(a*c*d)*(a+c)*(c-d)*(a+d)*y[2] ;
    const double C = -(a*b*d)*(a*b*d)*(a+b)*(a+d)*(-b+d)*y[3] ;
    const double D = (a*b*c)*(a*b*c)*(a+b)*(a+c)*(-b+c)*y[4] ;
    der[1] = ( A+B+C+D )/(a*-b*c*d*(a+b)*(a+c)*(-b+c)*(c-d)*(a+d)*(-b+d)) + \
      ( 1./a + 1./-b - 1./c - 1./d ) * y[1] ;
  }
  // quasi-symmetric solution, can be done in parallel
  size_t i ;
  for( i = 2 ; i < N-2 ; i++ ) {
    const double a = x[i] - x[i-2] ;
    const double b = x[i] - x[i-1] ;
    const double c = x[i+1] - x[i] ;
    const double d = x[i+2] - x[i] ;
    const double A = (b*c*d)*(b*c*d)*(b+c)*(c-d)*(b+d)*y[i-2] ;
    const double B = -(a*c*d)*(a*c*d)*(a+c)*(c-d)*(a+d)*y[i-1] ;
    const double C = -(a*b*d)*(a*b*d)*(a-b)*(a+d)*(b+d)*y[i+1] ;
    const double D = (a*b*c)*(a*b*c)*(a-b)*(a+c)*(b+c)*y[i+2] ;
    der[i] = ( A+B+C+D )/(a*b*c*d*(a-b)*(a+c)*(b+c)*(c-d)*(a+d)*(b+d)) + \
      ( 1./a + 1./b - 1./c - 1./d ) * y[i] ;
  }
  // -1 backward-looking O(h^4)
  {
    const double a = x[i] - x[i-3] ;
    const double b = x[i] - x[i-2] ;
    const double c = x[i] - x[i-1] ;
    const double d = x[i+1] - x[i] ;
    const double A = (b*c*d)*(b*c*d)*(b-c)*(-c-d)*(b+d)*y[i-3] ;
    const double B = -(a*c*d)*(a*c*d)*(a-c)*(-c-d)*(a+d)*y[i-2] ;
    const double C = -(a*b*d)*(a*b*d)*(a-b)*(a+d)*(b+d)*y[i-1] ;
    const double D = (a*b*c)*(a*b*c)*(a-b)*(a-c)*(b-c)*y[i+1] ;
    der[i] = ( A+B+C+D )/(a*b*-c*d*(a-b)*(a-c)*(b-c)*(-c-d)*(a+d)*(b+d)) + \
      ( 1./a + 1./b + 1./c - 1./d ) * y[i] ;
  }
  // -1 backward-looking O(h^4)
  {
    const double a = x[N-1] - x[N-5] ;
    const double b = x[N-1] - x[N-4] ;
    const double c = x[N-1] - x[N-3] ;
    const double d = x[N-1] - x[N-2] ;
    const double A = (b*c*d)*(b*c*d)*(b-c)*(-c+d)*(b-d)*y[N-5] ;
    const double B = -(a*c*d)*(a*c*d)*(a-c)*(-c+d)*(a-d)*y[N-4] ;
    const double C = -(a*b*d)*(a*b*d)*(a-b)*(a-d)*(b-d)*y[N-3] ;
    const double D = (a*b*c)*(a*b*c)*(a-b)*(a-c)*(b-c)*y[N-2] ;
    der[N-1] = ( A+B+C+D )/(a*b*-c*-d*(a-b)*(a-c)*(b-c)*(-c+d)*(a-d)*(b-d)) +\
      ( 1./a + 1./b + 1./c - 1./-d ) * y[N-1] ;
  }
  return ; // and that's all the derivatives
}

// cubic spline evaluator
double
cubic_eval( const double *__restrict x ,
	    const double *__restrict y ,
	    const double *__restrict der ,
	    const double mu ,
	    const int datalength )
{
  const size_t ref = find_idx( mu , x , datalength , 0 ) ;

  // newer, slightly higher order derivative
  const double diff = x[ref] - x[ref-1] ;
  const double yp   = der[ref-1] * diff ;
  const double yp_p = der[ref] * diff ;

  // coefficients for the poly y(x) = d + c*x + b*x^2 + a*x^3
  const double a = 2.0 * ( y[ref-1] - y[ref] ) + yp + yp_p ;
  const double b = 0.5 * ( -3.0 * a + ( yp_p - yp ) ) ;
  const double c = yp ;
  const double d = y[ref-1] ;

  // rescaled x-axis
  const double t = ( mu - x[ref-1] ) / diff ; // is the interpolating scale

  // Horner's cubic evaluation
  return d + t * ( c + t * ( b + t * a ) ) ;
}

// be careful, finds the minimum from the spline
double
cubic_min( const double *__restrict x ,
	   const double *__restrict y ,
	   const double *__restrict der ,
	   const size_t change_up )
{
  // if the cross over derivative is not bound we return 0
  if( change_up == 0 ) return 0.0 ;
 
  // finds the min as per evaluating the spline
  const double diff = x[change_up] - x[change_up-1] ;
  const double yp   = der[change_up-1] * diff ;
  const double yp_p = der[change_up] * diff ;

  const double a = 2.0 * ( y[change_up-1] - y[change_up] ) + yp + yp_p ;
  const double b = 0.5 * ( -3.0 * a + ( yp_p - yp ) ) ;
  const double c = yp ;

  // solve for y' = 0 :: 3ax^2 + 2bx + c = 0
  // should I be worried if a and b , or b and c or a and c are around 0? 
  // probably
  double ans = 0.0 ;
  if( fabs( a ) < 1E-15 ) {
    // edge case #1 a ~= 0, the rest are covered by below except the root
    ans = -c / ( 2.0 * b ) ;
  } else if( fabs( b ) < 1E-15 ) {
    // edge case #2 b ~= 0
    ans = sqrt( -c / ( 3.0 * a ) ) ;
  } else if( fabs( c ) < 1E-15 ) {
    // edge case #3 c ~= 0, non x=0 solution
    ans = -b / ( 1.5 * a ) ;
  } else {
    // should be a little more stable?
    const double root = fabs(b) * sqrt( 1.0 - 3. * a * c / ( b*b ) ) ;
    // positive root is the one we want
    ans = ( -b + root ) / ( 3.0 * a ) ;
  }

#ifdef verbose
  //Look at what we are solving ....
  fprintf( stdout , "[SPLINE] %f %f \n" , yp , yp_p ) ;
  fprintf( stdout , "[SPLINE] abc %e %e %e \n" , a , b , c ) ;
  fprintf( stdout , "[SPLINE] change up %zu \n" , change_up ) ;
  fprintf( stdout , "[SPLINE] root %e \n" , ans ) ;
#endif

  // and as we have solved for "t" in the note, we have to shift up
  return ans * diff + x[change_up-1] ;
}

// CALCULATES (one of) THE CUBE ROOT(s) //
inline static void 
cubert( double complex *__restrict res ,
	const double complex z )
{
  register const double real = cabs( z ) ;
  register const double angle = carg( z ) / 3.0 ; 
  *res = cbrt( real ) * ( cos( angle ) + I * sin( angle ) ) ; 
  return ;
}

// solve a quadratic ( ca x^2 + cb x + cc == 0 ) for x
static double
quad_solve( double complex *Z ,
	    const double ca ,
	    const double cb ,
	    const double cc ) 
{
  const double complex q = -0.5 * ( cb + sgn( cb ) * \
				    csqrt( cb*cb - 4.0 * ca*cc ) ) ;
  Z[ 0 ] = q / ca ;
  Z[ 1 ] = cc / q ;
  if( fabs( cimag( Z[0] ) ) < fabs( cimag( Z[0] ) ) ) {
    return creal( Z[0] ) ;
  } else {
    return creal( Z[1] ) ;
  }
}

// solve a cubic ( ca x^3 + cb x^2 + cc x + cd == eval ) for x 
double
cubic_solve( double complex Z[ 3 ] , 
	     const double ca ,
	     const double cb ,
	     const double cc , 
	     const double cd ,
	     const double eval )
{
  // if ca is super small this is a quadratic
  if( fabs( ca ) < PREC_TOL ) {
    return quad_solve( Z , cb , cc , cd - eval ) ; 
  }
  // if the constant term is zero the solution is quadratic too
  if( fabs( cd - eval ) < PREC_TOL ) {
    return quad_solve( Z , ca , cb , cc ) ; 
  }

  // divide through by a, x^3 + (b/a)x^2 + (c/a)x + ( d - eval )/ a
  const double a = ( cb / ca ) ;
  const double b = ( cc / ca ) ;
  const double c = ( ( cd - eval ) / ca ) ;

  const double Q = ( a * a - 3.0 * b ) / 9.0 ;
  const double R = ( a * ( 2.0 * a * a - 9.0 * b ) + \
		     27.0 * c ) / 54. ;

  // three real roots
  double theta ;
  if( R*R < Q*Q*Q ) {
    theta = acos( R / sqrt( Q*Q*Q ) ) ;
    Z[0] = -2.0 * sqrt( Q ) * cos( ( theta ) / 3.0 ) - a/3.0 ;
    Z[1] = -2.0 * sqrt( Q ) * cos( ( theta + TWOPI ) / 3.0 ) - a/3.0 ;
    Z[2] = -2.0 * sqrt( Q ) * cos( ( theta - TWOPI ) / 3.0 ) - a/3.0 ;
  } else {
    double complex res ;
    cubert( &res , R + csqrt( R*R - Q*Q*Q ) ) ;
    const double complex A = -res ;
    const double complex B = ( A != 0.0 ) ? Q/A : 0.0 ;
    Z[0] = ( A + B ) - a/3.0 ;
    Z[1] = -0.5 * ( A + B ) - a/3.0 - I * sqrt( 3 ) * 0.5 * ( A - B ) ;
    Z[2] = -0.5 * ( A + B ) - a/3.0 + I * sqrt( 3 ) * 0.5 * ( A - B ) ;
  }
  return 0.0 ;
}

// form our spline and solve at mu
double
solve_spline( const double *__restrict x ,
	      const double *__restrict y ,
	      const double *__restrict der ,
	      const double mu ,
	      const size_t change_up )
{
  if( change_up == 0 ) return -1.0 ;

  // newer, slightly higher order derivative
  const double diff = x[change_up] - x[change_up-1] ;
  const double yp   = der[change_up-1] * diff ;
  const double yp_p = der[change_up] * diff ;

  // coefficients for the poly y(x) = d + c*x + b*x^2 + a*x^3
  const double a = 2.0 * ( y[change_up-1] - y[change_up] ) + yp + yp_p ;
  const double b = 0.5 * ( -3.0 * a + ( yp_p - yp ) ) ;
  const double c = yp ;
  const double d = y[change_up-1] ;

  //return cubic_solve( a , b , c , d , mu ) ;//* diff + x[change_up-1] ;
  double complex Z[ 3 ] = { 0.0 , 0.0 , 0.0 } ;
  cubic_solve( Z , a , b , c , d , mu ) ;

  // should sanity-check the solution here too
  size_t roots = 0 ;
  for( roots = 0 ; roots < 3 ; roots++ ) {
    const double sol = creal( Z[ roots ] ) * diff + x[ change_up - 1 ] ;
    if( sol < x[ change_up ] && sol > x[ change_up - 1 ] ) {
      return sol ;
    }
  }
  return -1.0 ;
}

#undef TWOPI

#undef sgn
