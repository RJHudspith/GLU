/*
    Copyright 2013 Renwick James Hudspith

    This file (exactQ.c) is part of GLU.

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
   @file exactQ.c
   @brief logarithm routines and methods to compute the lie-fields of our link matrices
 */

#include "Mainfile.h"

#include "effs.h"     // computation of the f-constants
#include "solver.h"   // eigensolver using characteristic equation

// We have the bindings for Lapacke now
#ifdef HAVE_LAPACKE_H
  #include <lapacke.h> // we must be careful here
#else
// taylor approximation logs are pretty slow
  #include "taylor_logs.h"
#endif

// compute eigenvectors from s ,  put into v ,  need eigenvalues z//
static void 
vectors( GLU_complex *__restrict v ,
	 const GLU_complex *__restrict U ,
	 const double complex *__restrict z )
{
  GLU_complex S[ NCNC ] , B[ NCNC ] ; 
  size_t j , k , l ; 
  for( j = 0 ; j < NC ; j++ ) { // loop eigenvalues
    #if NC == 3
    S[0] = U[0] - (GLU_complex)z[j] ; S[1] = U[1] ; S[2] = U[2] ;
    S[3] = U[3] ; S[4] = U[4] - (GLU_complex)z[j] ; S[5] = U[5] ;
    S[6] = U[6] ; S[7] = U[7] ; S[8] = U[8] - (GLU_complex)z[j] ;
    #else
    size_t i ;
    for( i = 0 ; i < NCNC ; i++ ) {
      S[i] = ( i%( NC+1 ) == 0 ) ? U[i] - (GLU_complex)z[j] : U[i] ;
    }
    #endif
    // computes the adjunct of S //
    cofactor_transpose( B , S ) ; 
    //loop through B //
    for( k = 0 ; k < NC ; k++ ) {
      register GLU_real norm = 0. ;
      for( l = 0 ; l < NC ; l++ ) {
	const size_t place = k + NC * l ;
	norm += creal( B[place] ) * creal( B[place] ) +	\
	  cimag( B[place] ) * cimag( B[place] ) ;
      }
      if( norm > PREC_TOL ) {
	norm = 1.0 / sqrt( norm ) ;
	for( l = 0 ; l < NC ; l++ ) {
	  v[ j + NC * l  ] = B[ k + NC * l ] * norm ; 
	}
	break ;
      }
    }
  }
  return ;
}

// returns the ( A - A^{dag} ) * 0.5 of a matrix ( Antihermitian proj )
INLINE_VOID
AntiHermitian_proj( GLU_complex Q[ NCNC ] , 
		    const GLU_complex U[ NCNC ] )
{
#if NC == 3
  register GLU_real cimU0 = cimag( *( U + 0 ) ) ;
  register GLU_real cimU4 = cimag( *( U + 4 ) ) ;
  register GLU_real cimU8 = cimag( *( U + 8 ) ) ;
  *( Q + 0 ) = I * cimU0 ;
  *( Q + 1 ) = 0.5 * ( U[1] - conj( U[3] ) ) ;  
  *( Q + 2 ) = 0.5 * ( U[2] - conj( U[6] ) ) ; 
  *( Q + 3 ) = -conj( Q[1] ) ; 
  *( Q + 4 ) = I * cimU4 ;
  *( Q + 5 ) = 0.5 * ( U[5] - conj( U[7] ) ) ; 
  *( Q + 6 ) = -conj( Q[2] ) ;  
  *( Q + 7 ) = -conj( Q[5] ) ;  
  *( Q + 8 ) = I * cimU8 ;
#elif NC == 2
  register GLU_real cimU0 = cimag( *( U + 0 ) ) ;
  register GLU_real cimU3 = cimag( *( U + 3 ) ) ;
  *( Q + 0 ) = I * cimU0 ;
  *( Q + 1 ) = U[1] ; //0.5 * ( U[1] - conj( U[2] ) ) ;  
  *( Q + 2 ) = -conj( Q[1] )  ; 
  *( Q + 3 ) = I * cimU3 ;  
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Q[ j+i*NC ] = ( U[ j+i*NC ] - conj( U[ i+j*NC ] ) ) * 0.5 ;
    }
  }
#endif
  return ;
}

// Same as above but puts into the shortened hermitian form
INLINE_VOID
AntiHermitian_proj_short( GLU_complex Q[ HERMSIZE ] , 
			  const GLU_complex U[ NCNC ] ) 
{
#if NC == 3
  register GLU_real cimU0 = cimag( *( U + 0 ) ) ;
  register GLU_real cimU4 = cimag( *( U + 4 ) ) ;
  register GLU_real cimU8 = cimag( *( U + 8 ) ) ;
  *( Q + 0 ) = I * OneO3 * ( 2. * cimU0 - cimU4 - cimU8 ) ; 
  *( Q + 1 ) = 0.5 * ( U[1] - conj( U[3] ) ) ;  
  *( Q + 2 ) = 0.5 * ( U[2] - conj( U[6] ) ) ; 
  *( Q + 3 ) = I * OneO3 * ( 2. * cimU4 - cimU0 - cimU8 ) ;
  *( Q + 4 ) = 0.5 * ( U[5] - conj( U[7] ) ) ; 
#elif NC == 2
  register GLU_real cimU0 = cimag( *( U + 0 ) ) ;
  *( Q + 0 ) = I * cimU0 ;
  *( Q + 1 ) = U[1] ; //0.5 * ( U[1] - conj( U[2] ) ) ;  
#else
  size_t i , j , idx = 0 ;
  GLU_complex tr = 0. ;
  // compute the trace first
  for( i = 0 ; i < NC ; i++ ) {
    tr += cimag( U[ i*(NC+1) ] ) ;
  }
  tr /= (GLU_real)NC ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = i ; j < NC ; j++ ) { 
      Q[idx] = ( i != j ) ? ( U[ j + NC * i ] - conj( U[ i + NC * j ] ) ) * 0.5 : I * cimag( U[ j + NC * i ] ) - tr ;
      idx++ ;
    }
  }
#endif
  return ;
}

// compute the exact value of Q from principal log of e^( iQ ) from stefan
// durr's paper
void 
get_iQ( GLU_complex Q[ NCNC ] ,
	const GLU_complex U[ NCNC ] )
{
  size_t i , j ;
  double complex z[ NC ] ; 
  GLU_complex v[ NCNC ] , delta[ NCNC ] ; 
  Eigenvalues_suNC( z , U ) ; 
  vectors( v , U , z ) ; 
  // Finally compute V.Delta.V^{\dagger} //
  for( i = 0 ; i < NC ; i++ ) {
    const register double Z = carg( z[i] ) ;
    for( j = 0 ; j < NC ; j++ ) {
      delta[ j + i*NC ] = Z * conj( v[ i + j*NC ] ) ;
    }
  }
  multab( Q , v , delta ) ; 
  // test for problems with the trace being off by 2Pi
  #if NC == 3
  const GLU_real tr = creal( Q[0] ) + creal( Q[4] ) + creal( Q[8] ) ; 
  #elif NC == 2
  const GLU_real tr = creal( Q[0] ) + creal( Q[3] ) ; 
  #else
  const GLU_real tr = creal( trace ( Q ) ) ;
  #endif
  // recompute with a shifted evalue
  if( tr > PREC_TOL ) { // this one has tr = +2Pi
    register const double Z = carg( z[0] ) - TWOPI ;
    for( i = 0 ; i < NC ; i++ ) {
      delta[ i ] = Z * conj( v[ i*NC ] ) ;
    }
    multab( Q , v , delta ) ; 
  } else if ( creal( tr ) < -PREC_TOL  ) { // this one has tr = -2Pi
    register const double Z = carg( z[0] ) + TWOPI ;
    for( i = 0 ; i < NC ; i++ ) {
      delta[ i ] = Z * conj( v[ i*NC ] ) ;
    }
    multab( Q , v , delta ) ; 
  }
  return ;
}

// this code is called literally all the time, computes 
// Q in e^{iQ} up to roughly numerical precision and fairly quickly too
INLINE_VOID
exact_log_slow( GLU_complex Q[ NCNC ] , 
		const GLU_complex U[ NCNC ] ) 
{
#if NC == 3   
  double complex f[ NC ] , z[ NC ] ;
  Eigenvalues_suNC( z , U ) ; 
  z[0] = carg( z[0] ) ; 
  z[1] = carg( z[1] ) ;
  // This is now "numerically stable", just as the MP prescription
  // logically the log of these gives the hermitian z's
  z[2] = -( z[0] - z[1] ) ; 
  if( unlikely( z[0] == z[1] ) && unlikely( z[1] == z[2] ) ) {
    return ;
  }
  f_hermitian_log_suNC( f , z ) ;

  // I think we can do better
  const double complex con = conj( f[2] ) ;
  const double mod = 1.0 / ( cimag( f[1] ) * creal( f[2] ) -
			     creal( f[1] ) * cimag( f[2] ) ) ;
  const double imf = ( cimag( f[0] ) * creal( f[2] ) -
		       creal( f[0] ) * cimag( f[2] ) ) * mod ;
  const double complex trce = mod * OneOI2 ; 

  // complete the resulting hermitian matrix, 
  // "forcing" it slightly to be hermitian
  // to alleviate any resulting round-off errors.
  *( Q + 0 ) = trce * ( con * U[0] - f[2] * conj( U[0] ) ) - imf ; 
  *( Q + 1 ) = trce * ( con * U[1] - f[2] * conj( U[3] ) ) ; 
  *( Q + 2 ) = trce * ( con * U[2] - f[2] * conj( U[6] ) ) ; 
  *( Q + 3 ) = conj( Q[1] ) ;  
  *( Q + 4 ) = trce * ( con * U[4] - f[2] * conj( U[4] ) ) - imf ;  
  *( Q + 5 ) = trce * ( con * U[5] - f[2] * conj( U[7] ) ) ; 
  *( Q + 6 ) = conj( Q[2] ) ; 
  *( Q + 7 ) = conj( Q[5] ) ;  
  *( Q + 8 ) = -( Q[0] + Q [4] ) ; 
#elif NC == 2
  // inlined the two necessary routines ...
  register const double reU0 = (double)creal( U[0] ) ;
  const double complex z = reU0 + I * sqrt ( 1.0 - reU0 * reU0 ) ;
  const double herm_z = carg( z ) ;
  const double f0 = creal( z ) ; //cos( herm_z ) ;
  const double f1 = ( fabs(herm_z) < SINTOLSU2 ) ? ( 1. - herm_z / 6. * ( 1. - herm_z / 20. * ( 1. - herm_z / 42. ) ) ) : cimag( z ) / herm_z ;

  // use Cayley-Hamilton again A = ( U - f0 ) / f1 //
  // for su2 this is nice as they both tend to 1
  // this of course pushes the value of A -> 0
  register const double complex oneOf1 = -I / ( f1 ) ;
  *( Q + 0 ) = ( U[0] - f0 ) * oneOf1 ;
  *( Q + 1 ) = ( U[1] ) * oneOf1 ;
  *( Q + 2 ) = conj( Q[1] ) ;
  *( Q + 3 ) = -Q[0] ;
#else
  #ifdef HAVE_LAPACKE_H
  GLU_complex b[ NCNC ] , evalues[ NC ] , *vl , vr[ NCNC ] ;
  memcpy( b , U , NCNC * sizeof( GLU_complex ) ) ;
  const size_t n = NC , lda = NC  , ldvl = NC , ldvr = NC ;
  int info = LAPACKE_zgeev( LAPACK_ROW_MAJOR , 'N' , 'V' ,
			    n , b , lda , evalues , 
			    vl, ldvl, 
			    vr, ldvr ) ;
  // something broke here
  if( info != 0 ) { 
    fprintf( stderr , "[EXACTQ] info :: %d \n" , info ) ;
    size_t i ;
    for( i = 0 ; i < NC ; i++ ) { printcomplex( evalues[i] ) ; }
  }
  // complete with V.D.V^{\dagger}
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    register const GLU_real zi = carg( evalues[i] ) ;
    for( j = 0 ; j < NC ; j++ ) {
      b[ j + i * NC ] = conj( vr[ i + j*NC ] ) * zi ;
    }   
  }
  multab( Q , vr , b ) ;
  #else
  // get desperate and call the series log
  brute_force_log( Q , U , NC ) ;
#endif
#endif
  return ;
}

// performs that of above but for the shortened Hermitian matrix form
INLINE_VOID
exact_log_slow_short( GLU_complex Q[ HERMSIZE ] , 
		      const GLU_complex U[ NCNC ] )
{
#if NC == 3 
  double complex f[NC] , z[NC] ;
  Eigenvalues_suNC( z , U ) ; 
  z[0] = carg( z[0] ) ; 
  z[1] = carg( z[1] ) ;
  // This is now "numerically stable", just as the MP prescription
  // logically the log of these gives the hermitian z's
  z[2] = -( z[0] - z[1] ) ; 
  register const double a = 0.5 * creal( z[0] ) ;
  const double u = a ;
  register const double w = creal( z[1] ) + a ;
  const double cw = cos( w ) ;
  register const double cu = cos( u ) ;
  register const double su = sin( u ) ; //sin( u ) ;
  const double complex one1 = cu - I * su ;
  double complex two1 = conj( one1 ) ; //cu + I * su ;
  two1 *= two1 ;
  const double uu = u * u ;
  const double ww = w * w ; 
  // control for degeneracies, if we are degenerate the solution is basically
  // the identity matrix, is definitely an ad hoc solution
  const double denom = 1.0 / ( 9. * uu - ww ) ;
  //double E0 = 0. ; 
  // we only allow the taylor expansion at a very low acc to allow the 3D
  // hyl smearing to not get stuck, this is similar to the problem I saw with
  // the su2 log and stout expansions
  const double E0 = fabs( w ) < STOL ? 1 - ww / 6. * ( 1 - ww / 20. * ( 1 - ww / 42. ) ) : sin( w ) / w ; 
  //const register double twou = 2. * u ;
  f[0] = ( uu - ww ) * two1 + one1 * ( 8 * uu * cw + I * 2 * u * ( 3 * uu + ww ) * E0 ) ; 
  f[1] = 2 * u * two1 - one1 * ( 2 * u * cw - I * ( 3 * uu - ww ) * E0 ) ; 
  f[2] = two1 - one1 * ( cw + 3 * I * u * E0 ) ; 
  f[0] *= denom ; 
  f[1] *= denom ; 
  f[2] *= denom ; 
  const double complex con = conj( f[2] ) ;
  const double mod = 1.0 / cimag( f[1] * con ) ;
  const double imf = cimag( f[0] * con ) * mod ; 
  const double complex trce = mod * OneOI2 ; 
  // complete the resulting hermitian matrix, 
  // "forcing" it slightly to be hermitian
  // to alleviate any resulting round-off errors.
  *( Q + 0 ) = trce * ( con * U[0] - f[2] * conj( U[0] ) ) - imf ; 
  *( Q + 1 ) = trce * ( con * U[1] - f[2] * conj( U[3] ) ) ; 
  *( Q + 2 ) = trce * ( con * U[2] - f[2] * conj( U[6] ) ) ; 
  *( Q + 3 ) = trce * ( con * U[4] - f[2] * conj( U[4] ) ) - imf ;  
  *( Q + 4 ) = trce * ( con * U[5] - f[2] * conj( U[7] ) ) ; 
#elif NC == 2
  register const double reU0 = (double)creal( U[0] ) ;
  const double complex z = reU0 + I * sqrt ( 1.0 - reU0 * reU0 ) ;
  const double herm_z = carg( z ) ;
  const double f0 = reU0 ; //cos( herm_z ) ;
  const double f1 = ( fabs( herm_z ) < SINTOLSU2 ) ? ( 1 - herm_z / 6. * ( 1 - herm_z / 20. * ( 1 - herm_z / 42. ) ) ) : cimag( z ) / herm_z ;
  register const double complex oneOf1 = -I / ( f1 ) ;
  // use Cayley-Hamilton again A = ( U - f0 ) / f1 //
  // for su2 this is nice as they both tend to 1
  // this of course pushes the value of A -> 0
  *( Q + 0 ) = ( U[0] - f0 ) * oneOf1 ;
  *( Q + 1 ) = ( U[1] ) * oneOf1 ;
#else
  GLU_complex QP[ NCNC ] ;
  exact_log_slow( QP , U ) ;
  pack_hermitian( Q , QP ) ;
#endif
  return ;
}

// the hermitian projection of the logarithm
INLINE_VOID
Hermitian_proj( GLU_complex Q[ NCNC ] , 
		const GLU_complex U[ NCNC ] )
{
#if NC == 3
  register const GLU_real cimU0 = cimag( *( U + 0 ) ) ;
  register const GLU_real cimU4 = cimag( *( U + 4 ) ) ;
  register const GLU_real cimU8 = cimag( *( U + 8 ) ) ;
  *( Q + 0 ) = OneO3 * ( 2 * cimU0 - cimU4 - cimU8 ) ; 
  *( Q + 1 ) = OneOI2 * ( U[1] - conj( U[3] ) ) ;  
  *( Q + 2 ) = OneOI2 * ( U[2] - conj( U[6] ) ) ; 
  *( Q + 3 ) = conj( Q[1] ) ; 
  *( Q + 4 ) = OneO3 * ( 2 * cimU4 - cimU0 - cimU8 ) ;
  *( Q + 5 ) = OneOI2 * ( U[5] - conj( U[7] ) ) ; 
  *( Q + 6 ) = conj( Q[2] ) ;  
  *( Q + 7 ) = conj( Q[5] ) ;  
  *( Q + 8 ) = -creal( Q[0] ) - creal( Q[4] ) ; 
#elif NC == 2
  *( Q + 0 ) = cimag( U[0] ) ;
  *( Q + 1 ) = -I * U[1] ; //OneOI2 * ( U[1] - conj( U[2] ) ) ;  
  *( Q + 2 ) = conj( Q[1] )  ; 
  *( Q + 3 ) = -Q[0] ;  
#else
  int i , j ;
  register double tr = 0.0 ;
  for( i = 0 ; i < NC ; i++ ) {
    Q[ i*(NC+1) ] = cimag( U[ i*(NC+1) ] ) ;
    tr += creal( Q[ i*(NC+1) ] ) ; 
    for( j = i+1 ; j < NC ; j++ ) {
      Q[ j+i*NC ] = ( U[ j+i*NC ] - conj( U[ i+j*NC ] ) ) * OneOI2 ;
      Q[ i+j*NC ] = conj( Q[ j+i*NC ] ) ;
    }
  }
  tr /= NC ;
  for( i = 0 ; i < NC ; i++ ) {
    Q[ i*(NC+1) ] -= tr ;
  }
#endif
  return ;
}

// Same as above but puts into the shortened hermitian form
INLINE_VOID
Hermitian_proj_short( GLU_complex Q[ HERMSIZE ] , 
		      const GLU_complex U[ NCNC ] )
{
#if NC == 3
  register const GLU_real cimU0 = cimag( *( U + 0 ) ) ;
  register const GLU_real cimU4 = cimag( *( U + 4 ) ) ;
  register const GLU_real cimU8 = cimag( *( U + 8 ) ) ;
  *( Q + 0 ) = OneO3 * ( 2. * cimU0 - cimU4 - cimU8 ) ; 
  *( Q + 1 ) = OneOI2 * ( U[1] - conj( U[3] ) ) ;  
  *( Q + 2 ) = OneOI2 * ( U[2] - conj( U[6] ) ) ; 
  *( Q + 3 ) = OneO3 * ( 2. * cimU4 - cimU0 - cimU8 ) ;
  *( Q + 4 ) = OneOI2 * ( U[5] - conj( U[7] ) ) ; 
#elif NC == 2
  *( Q + 0 ) = cimag( U[0] ) ;
  *( Q + 1 ) = -I * U[1] ; //OneOI2 * ( U[1] - conj( U[2] ) ) ;  
#else
  size_t i , j , idx = 0 ;
  register double tr = 0. ;
  // compute the trace first
  for( i = 0 ; i < NC ; i++ ) {
    tr += cimag( U[ i*(NC+1) ] ) ;
  }
  tr /= (double)NC ;
  // fill up "Q"
  for( i = 0 ; i < NC-1 ; i++ ) {
    Q[idx] = cimag( U[ i*(NC+1) ] ) - tr ;
    idx++ ;
    for( j = i+1 ; j < NC ; j++ ) { 
      Q[idx] = ( U[ j + NC * i ] - conj( U[ i + NC * j ] ) ) * OneOI2 ;
      idx++ ;
    }
  }
#endif
  return ;
}

// returns the traceless ( A - A^{dag} ) * 0.5 of a matrix ( Antihermitian proj )
INLINE_VOID
trf_AntiHermitian_proj( GLU_complex Q[ NCNC ] , 
			const GLU_complex U[ NCNC ] )
{
#if NC == 3
  register GLU_real cimU0 = cimag( *( U + 0 ) ) ;
  register GLU_real cimU4 = cimag( *( U + 4 ) ) ;
  register GLU_real cimU8 = cimag( *( U + 8 ) ) ;
  *( Q + 0 ) = I * OneO3 * ( 2. * cimU0 - cimU4 - cimU8 ) ;
  *( Q + 1 ) = 0.5 * ( U[1] - conj( U[3] ) ) ;  
  *( Q + 2 ) = 0.5 * ( U[2] - conj( U[6] ) ) ; 
  *( Q + 3 ) = -conj( Q[1] ) ; 
  *( Q + 4 ) = I * OneO3 * ( 2. * cimU4 - cimU0 - cimU8 ) ;
  *( Q + 5 ) = 0.5 * ( U[5] - conj( U[7] ) ) ; 
  *( Q + 6 ) = -conj( Q[2] ) ;  
  *( Q + 7 ) = -conj( Q[5] ) ;  
  *( Q + 8 ) = -Q[0] - Q[4] ;
#elif NC == 2
  register GLU_real cimU0 = cimag( *( U + 0 ) ) ;
  register GLU_real cimU3 = cimag( *( U + 3 ) ) ;
  *( Q + 0 ) = I * cimU0 ;
  *( Q + 1 ) = U[1] ; //0.5 * ( U[1] - conj( U[2] ) ) ;  
  *( Q + 2 ) = -conj( Q[1] )  ; 
  *( Q + 3 ) = I * cimU3 ;  
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Q[ j+i*NC ] = ( U[ j+i*NC ] - conj( U[ i+j*NC ] ) ) * 0.5 ;
    }
  }
  GLU_complex tr = trace( Q ) / (GLU_real)NC ;
  for( i = 0 ; i < NC ; i++ ) {
    Q[ i*(NC+1) ] -= tr ;
  }
#endif
  return ;
}
