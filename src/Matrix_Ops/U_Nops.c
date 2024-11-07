/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (U_Nops.c) is part of GLU.

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
   @file U_Nops.c
   @brief matrix operations
 */
#include "Mainfile.h"

#include "LU.h" // determinant

// add constant diagonal matrix to a
inline void
add_constant( GLU_complex a[ NCNC ] , 
	      const GLU_complex c ) 
{
#if NC == 3
  a[0] += c ;
  a[4] += c ;
  a[8] += c ;
#elif NC == 2
  a[0] += c ;
  a[3] += c ;
#else
  size_t i ;
  for( i = 0 ; i < NC ; i++ ) {
    a[ i*(NC+1) ] += c ;
  }
#endif
  return ;
}

// (atomic) a += b //
void
a_plus_b( GLU_complex a[ NCNC ] , 
	  const GLU_complex b[ NCNC ] )
{
#if NC == 3
  a[0] += b[0] ;  a[1] += b[1] ;  a[2] += b[2] ; 
  a[3] += b[3] ;  a[4] += b[4] ;  a[5] += b[5] ; 
  a[6] += b[6] ;  a[7] += b[7] ;  a[8] += b[8] ; 
#elif NC == 2
  a[0] += b[0] ;  a[1] += b[1] ;  
  a[2] += b[2] ;  a[3] += b[3] ; 
#else
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[i] += b[i] ;
  }
#endif
  return ;
}

// a += S.b where S is complex //
void
a_plus_CSxb( GLU_complex a[ NCNC ] , 
	     const GLU_complex b[ NCNC ] ,
	     const GLU_complex S )
{
#if NC == 3
  a[0] += S * b[0] ;  a[1] += S * b[1] ;  a[2] += S * b[2] ; 
  a[3] += S * b[3] ;  a[4] += S * b[4] ;  a[5] += S * b[5] ; 
  a[6] += S * b[6] ;  a[7] += S * b[7] ;  a[8] += S * b[8] ; 
#elif NC == 2
  a[0] += S * b[0] ;  a[1] += S * b[1] ;  
  a[2] += S * b[2] ;  a[3] += S * b[3] ; 
#else
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[i] += S * b[i] ;
  }
#endif
  return ;
}

// (atomic) a += S.b //
void
a_plus_Sxb( GLU_complex a[ NCNC ] , 
	    const GLU_complex b[ NCNC ] ,
	    const GLU_real S )
{
#if NC == 3
  a[0] += S * b[0] ;  a[1] += S * b[1] ;  a[2] += S * b[2] ; 
  a[3] += S * b[3] ;  a[4] += S * b[4] ;  a[5] += S * b[5] ; 
  a[6] += S * b[6] ;  a[7] += S * b[7] ;  a[8] += S * b[8] ; 
#elif NC == 2
  a[0] += S * b[0] ;  a[1] += S * b[1] ;  
  a[2] += S * b[2] ;  a[3] += S * b[3] ; 
#else
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[i] += S * b[i] ;
  }
#endif
  return ;
}

// atomically add short matrices A += S * ( B - C )
void
a_plus_Sxbminc_short( GLU_complex a[ HERMSIZE ] ,
		      const GLU_real S ,
		      const GLU_complex B[ HERMSIZE ] ,
		      const GLU_complex C[ HERMSIZE ] )
{
#if NC == 3
  a[0] += S * ( B[ 0 ] - C[ 0 ] ) ;
  a[1] += S * ( B[ 1 ] - C[ 1 ] ) ;
  a[2] += S * ( B[ 2 ] - C[ 2 ] ) ;
  a[3] += S * ( B[ 3 ] - C[ 3 ] ) ;
  a[4] += S * ( B[ 4 ] - C[ 4 ] ) ;
#elif NC == 2 
  a[0] += S * ( B[ 0 ] - C[ 0 ] ) ;
  a[1] += S * ( B[ 1 ] - C[ 1 ] ) ;
#else
  size_t mu ;
  for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
    a[ mu ] += S * ( B[ mu ] - C[ mu ] ) ;
  }
#endif
  return ;
}

// matrix a = b - c //
void
b_min_c( GLU_complex a[ NCNC ] , 
	 const GLU_complex b[ NCNC ] ,
	 const GLU_complex c[ NCNC ] )
{
#if NC == 3
  a[0] = b[0] - c[0] ; a[1] = b[1] - c[1] ; a[2] = b[2] - c[2] ; 
  a[3] = b[3] - c[3] ; a[4] = b[4] - c[4] ; a[5] = b[5] - c[5] ; 
  a[6] = b[6] - c[6] ; a[7] = b[7] - c[7] ; a[8] = b[8] - c[8] ; 
#elif NC == 2
  a[0] = b[0] - c[0] ;  a[1] = b[1] - c[1] ;  
  a[2] = b[2] - c[2] ;  a[3] = b[3] - c[3] ;  
#else
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[i] = b[i] - c[i] ;
  }
#endif
  return ;
}

// returns the determinant
GLU_complex
cofactor_transpose( GLU_complex a[ NCNC ] ,
		    const GLU_complex b[ NCNC ] )
{
#if NC == 3
  a[0] = b[4] * b[8] - b[5] * b[7] ; 
  a[1] = -( b[1] * b[8] - b[2] * b[7] ) ; 
  a[2] = b[1] * b[5] - b[2] * b[4] ;
  a[3] = -( b[3] * b[8] - b[5] * b[6] ) ;
  a[4] = b[0] * b[8] - b[2] * b[6] ; 
  a[5] = -( b[0] * b[5] - b[2] * b[3] ) ; 
  a[6] = b[3] * b[7] - b[4] * b[6] ; 
  a[7] = -( b[0] * b[7] - b[1] * b[6] ) ; 
  a[8] = b[0] * b[4] - b[1] * b[3] ; 
  return b[0] * a[0] + b[1] * a[3] + b[2] * a[6] ;
#elif NC == 2
  a[0] = b[3] ;
  a[1] = -b[1] ;
  a[2] = -b[2] ;
  a[3] = b[0] ;
  return b[0] * a[0] + b[1] * a[2] ;
#else
  // compute the cofactor matrix
  GLU_complex array[ ( NC - 1 ) * ( NC - 1 ) ]  GLUalign , 
    temp[ NCNC ] GLUalign ;
  size_t i , j ;
  for( i = 0 ; i < NCNC ; i++ ) { // our bona-fide minor index
    size_t idx = 0 ;
    for( j = 0 ; j < ( NC * NC ) ; j++ ) {
      if( ( j%NC != i%NC ) && ( j/NC != i/NC ) ) { // remove columns and rows
	// pack array
	array[idx] = b[j] ;
	idx ++ ;
      } 
    }
    // compute the determinant
    // factor is (-1)^{row_idx+c_idx}
    const GLU_real mulfact = ( i/NC + i%NC )%2 ? -1.0 : 1.0 ; 
    temp[i] = mulfact * (GLU_complex)LU_det( NC - 1 , array ) ;
  }
  // compute the determinant
  GLU_complex det = 0.0 ;
  for( i = 0 ; i < NC ; i++ ) {
    det += b[ i ] * temp[ i ] ; 
  }
  transpose( a , temp ) ;
  return det ;
#endif
}

//computes the transpose and conjugate of the matrix a put into b //
void
dagger( GLU_complex b[ NCNC ] , 
	const GLU_complex a[ NCNC ] )
{
#if NC == 3
  b[0] = conj( a[0] ) ; b[1] = conj( a[3] ) ; b[2] = conj( a[6] ) ; 
  b[3] = conj( a[1] ) ; b[4] = conj( a[4] ) ; b[5] = conj( a[7] ) ; 
  b[6] = conj( a[2] ) ; b[7] = conj( a[5] ) ; b[8] = conj( a[8] ) ; 
#elif NC == 2
  b[0] = conj( a[0] )  ;  b[1] = conj( a[2] ) ; 
  b[2] = conj( a[1] )  ;  b[3] = conj( a[3] ) ;  
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      b[ j + i * NC ] = conj( a[ i + j * NC ] ) ;
    }
  }
#endif
  return ;
}

// Determinant //
GLU_complex
det( const GLU_complex U[ NCNC ] ) 
{
#if NC == 3
  return U[0] * ( U[4] * U[8]-U[5] * U[7] ) -	\
    U[1] * ( U[3] * U[8] - U[5] * U[6] ) +	\
    U[2] * ( U[3] * U[7] - U[4] * U[6] ) ; 
#elif NC == 2
  return U[0] * U[3] - U[1] * U[2] ; 
#else
  return (GLU_complex)LU_det( NC , U ) ;
#endif
}

// diagonal matrix of some value "c" //
void
diag( GLU_complex M[ NCNC ] ,
      const GLU_complex c )
{
#if NC == 3
  M[0] = c  ; M[1] = 0. ; M[2] = 0. ; 
  M[3] = 0. ; M[4] = c  ; M[5] = 0. ; 
  M[6] = 0. ; M[7] = 0. ; M[8] = c  ; 
#elif NC == 2
  M[0] = c  ; M[1] = 0. ; 
  M[2] = 0. ; M[3] = c  ;
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      M[ j + NC * i ] = ( i != j ) ? 0.0 : c ;
    }
  }
#endif
  return ;
}

// matrix with an NC-vector on the diagonal //
void
diag_vect( GLU_complex M[ NCNC ] ,
	   const GLU_complex c[ NC ] )
{
#if NC == 3
  M[0] = c[0] ; M[1] = 0.   ; M[2] = 0.   ; 
  M[3] = 0.   ; M[4] = c[1] ; M[5] = 0.   ; 
  M[6] = 0.   ; M[7] = 0.   ; M[8] = c[2] ; 
#elif NC == 2
  M[0] = c[ 0 ] ; M[1] = 0.     ; 
  M[2] = 0.     ; M[3] = c[ 1 ] ;
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = i+1 ; j < NC ; j++ ) {
      M[ j + i * NC ] = M[ i + j * NC ] = 0.0 ;
    }
    M[ i * ( NC + 1 ) ] = c[ i ] ;
  }
#endif
  return ;
}

//equivalence set a equal to b//
void
equiv( GLU_complex a[ NCNC ] , 
       const GLU_complex b[ NCNC ] )
{
#if NC == 3
  a[0] = b[0] ;  a[1] = b[1] ;  a[2] = b[2] ; 
  a[3] = b[3] ;  a[4] = b[4] ;  a[5] = b[5] ; 
  a[6] = b[6] ;  a[7] = b[7] ;  a[8] = b[8] ; 
#elif NC == 2
  a[0] = b[0] ;  a[1] = b[1] ;  
  a[2] = b[2] ;  a[3] = b[3] ;  
#else
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = b[ i ] ;
  }
#endif
  return ;
}

// NCxNC identity matrix in complex //
void
identity( GLU_complex ident[ NCNC ] )
{
#if NC == 3
  ident[0] = 1. ; ident[1] = 0. ; ident[2] = 0. ; 
  ident[3] = 0. ; ident[4] = 1. ; ident[5] = 0. ; 
  ident[6] = 0. ; ident[7] = 0. ; ident[8] = 1. ; 
#elif NC == 2
  ident[0] = 1.0  ;  ident[1] = 0. ; 
  ident[2] = 0.0  ;  ident[3] = 1.0 ; 
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ){
    for( j = 0 ; j < NC ; j++ ) {
      ident[ j + NC * i ] = ( i != j ) ? 0.0 : 1.0 ;
    }
  }
#endif
  return ;
}

// check if it IS UNITARY //
GLU_bool
is_unitary( const GLU_complex U[ NCNC ] ) 
{
  GLU_complex temp[ NCNC ] GLUalign ; 
  GLU_real vv = 0. ; 
  GLU_bool unitary = GLU_TRUE ;
  size_t i , j ; 

  // check U.U^{dag} = I
  multab_dag( temp , U , U ) ; 

  //first do U.U^{dagger} and check it is unitary
  for( i = 0 ; i < NCNC ; i++ ) {
    const GLU_real check = !( i%(NC+1) ) ? cabs( temp[i] ) - 1.0 :\
      cabs( temp[i] ) ; 
    if( fabs( check ) > PREC_TOL ) {
      unitary = GLU_FALSE ;
      fprintf( stderr , "* flag seen * element -> %zu :: %1.8e\n" , 
	       i , check ) ;
      fprintf( stderr , "%1.15f %1.15f \n" , 
	       creal( temp[i] ) , cimag( temp[i] ) ) ;
      goto end ;
    }
    // same as isnan()
    if( U[i] != U[i] ) {
      unitary = GLU_FALSE ;
      fprintf( stderr , "* We have a NaN here * element %zu %f %f\n" , 
	       i , creal( U[i] ) , cimag( U[i] ) ) ;
      goto end ;
    }
  }

  //also check orthogonality of columns and rows
  //columns
  for( j = 0 ; j < NC ; j++ ) {
    for( i = 0 ; i < NC ; i++ ) {
      vv += creal( U[ NC * i ] * conj( U[ NC * i + j ] ) ) ; 
    }
    if( fabs( vv - 1.0 )/NC > PREC_TOL ) {
      unitary = GLU_FALSE ;
      fprintf( stderr , "[column %zu] not orthogonal!!!\n"
	       "MUST REUNITARIZE [reunit( Z , U )] vv %1.8f" , j , vv ) ; 
      goto end ;
    }
  }

  vv = 0 ; 
  //check if rows are orthogonal
  for( j = 0 ; j < NC ; j++ ) {
    for( i = 0 ; i < NC ; i++ ) {
      vv += creal( U[ i ] * conj( U[ i + NC * j ] ) ) ; 
    }
    if( fabs( vv - 1.0 )/NC > PREC_TOL ) {
      unitary = GLU_TRUE ;
      fprintf( stderr , "[row %zu] not orthogonal!!!\nMUST REUNITARIZE"
	       " [reunit( Z , U )] vv %1.8f" , j , vv ) ; 
      goto end ;
    }
  }

  // if we see a problem we print out what should be the identity


 end :
  if( unitary == GLU_FALSE ) { write_matrix( temp ) ; }

  return unitary ;
}

// matrix times a vector ( vector )vect=( matrix )S.( vector )v //
void
mat_mult_vec( GLU_complex vect[ NC ] , 
	      const GLU_complex S[ NCNC ] , 
	      const GLU_complex v[ NC ] ) 
{
#if NC == 3
  vect[0] = S[0] * v[0] + S[1] * v[1] + S[2] * v[2] ; 
  vect[1] = S[3] * v[0] + S[4] * v[1] + S[6] * v[2] ; 
  vect[2] = S[6] * v[0] + S[7] * v[1] + S[8] * v[2] ; 
#elif NC == 2
  vect[0] = S[ 0 ] * v[ 0 ] + S[ 1 ] * v[ 1 ] ;
  vect[1] = S[ 2 ] * v[ 0 ] + S[ 3 ] * v[ 1 ] ;
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    vect[ i ] = 0. ;
    for( j = 0 ; j < NC ; j++ ) {
      vect[i] += S[ j + i * NC ] * v[ j ] ;  
    }
  }
#endif
  return ;
}

// ( M )atrix multiplied by a constant //
void
M_times_c( GLU_complex M[ NCNC ] , 
	   const GLU_complex c ) 
{
#if NC == 3
  M[0] *= c ; M[1] *= c ; M[2] *= c ; 
  M[3] *= c ; M[4] *= c ; M[5] *= c ; 
  M[6] *= c ; M[7] *= c ; M[8] *= c ; 
#elif NC == 2
  M[0] *= c ; M[1] *= c ; 
  M[2] *= c ; M[3] *= c ; 
#else
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    M[ i ] *= c ;
  }
#endif
  return ;
}

// nodes for AntiHermitian_proj linked lists
struct node {
  GLU_bool squareable ;
  struct node *next ;
} ;

// matrix power routine using a simple linked list structure
void
matrix_power( GLU_complex a[ NCNC ] , 
	      const GLU_complex b[ NCNC ] , 
	      const int n )
{
  if( n == 0 ) { identity( a ) ; return ; } 
  else if( n == 1 ) { equiv( a , b ) ; return ; } 
  else if( n == 2 ) { multab( a , b , b ) ; return ; } 
  else {
    // generate our linked list
    struct node *head = NULL , *curr ;
    size_t nn = n , length = 0 , i ;
    // compute the list, and its length
    while( nn > 2 ) {
      curr = (struct node*)malloc( sizeof( struct node ) ) ;
      if( nn%2 == 0 ) {
	nn = nn >> 1 ;
	curr -> squareable = GLU_TRUE ;
      } else {
	nn -- ;
	curr -> squareable = GLU_FALSE ;
      }
      length ++ ;
      curr -> next = head ;
      head = curr ;
    }
    curr = head ;
    // go back through the list performing squarings if possible
    GLU_complex tmp[ NCNC ] GLUalign ;
    multab( a , b , b ) ; 
    for( i = 0 ; i < length ; i++ ) {
      if( curr -> squareable != GLU_FALSE ) {
	multab( tmp , a , a ) ; 
      } else {
	multab( tmp , a , b ) ;
      }
      equiv( a , tmp ) ;
      curr = curr -> next ;
    }
    // clean up the list, removing all the allocs
    while( ( curr = head ) != NULL ) {
      head = head -> next ;
      free( curr ) ;
    }
  }
  return ;
}

// computes a[NC]*b[NC]^{dagger} puts into Q[NCNC] ( OuterProduct ) //
void
outerproduct( GLU_complex Q[ NCNC ] , 
	      const GLU_complex a[ NC ] , 
	      const GLU_complex b[ NC ] )
{
#if NC == 3
  Q[0] = a[0] * conj( b[0] ) ; Q[1] = a[0] * conj( b[1] ) ; Q[2] = a[0] * conj( b[2] ) ; 
  Q[3] = a[1] * conj( b[0] ) ; Q[4] = a[1] * conj( b[1] ) ; Q[5] = a[1] * conj( b[2] ) ; 
  Q[6] = a[2] * conj( b[0] ) ; Q[7] = a[2] * conj( b[1] ) ; Q[8] = a[2] * conj( b[2] ) ; 
#elif NC == 2
  Q[0] = a[0] * conj( b[0] ) ; Q[1] = a[0] * conj( b[1] ) ; 
  Q[2] = a[1] * conj( b[0] ) ; Q[3] = a[1] * conj( b[1] ) ; 
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Q[ j + i * NC ] = a[ i ] * conj( b[ j ] ) ;
    }
  }
#endif
  return ;
}

// pack a hermitian matrix into the "HERMSIZE" packing
void
pack_hermitian( GLU_complex a[ HERMSIZE ] ,
		const GLU_complex b[ NCNC ] )
{
  // pack into hermsize
  size_t i , j , idx = 0 ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = i ; j < NC ; j++ ) {
      a[ idx ] = b[ j + i*NC ] ;
      idx++ ;
    }
  }
  return ;
}

// print complex function
void
printcomplex( const GLU_complex a ) 
{
  fprintf( stdout , " ( %1.7e  ,  %1.7e )\n" , creal( a ) , cimag( a ) ) ; 
  return ;
}

#if NC<4
// rebuild the 8 parameter configuration
void 
rebuild( GLU_complex a[ NCNC ] , 
	 const GLU_real *__restrict b )
{
#if NC == 3
  const double tmp = ( b[0] * b[0] + b[1] * b[1] + b[2] * b[2] + b[3] * b[3] ) ;
  #ifdef SINGLE_PREC
  const GLU_bool test = ( fabs( tmp ) < FLT_MIN ) ? GLU_TRUE : GLU_FALSE ;
  #else
  const GLU_bool test = ( fabs( tmp ) < DBL_MIN ) ? GLU_TRUE : GLU_FALSE ;
  #endif
  if( test != GLU_FALSE ) {
    identity( a ) ;
  } else {
    const double NORM = 1.0 / ( tmp ) ; 
    // be super careful about rounding precisions ...
    register const double reg1 = 1.0 - ( b[0] * b[0] + b[1] * b[1] ) - ( b[2] * b[2] + b[3] * b[3] ) ;
    // create the variable a1 ; 
    const double complex a1 = reg1 < 0.0 ? 0.0 : sqrt( reg1 ) * ( cos( b[6] ) + I * sin( b[6] ) ) ; 
    register const double reg2 = 1.0 - ( creal(a1)*creal(a1) + cimag(a1)*cimag(a1) ) - ( b[4] * b[4] + b[5] * b[5] ) ;
    const double complex c = reg2 < 0.0 ? 0.0 : sqrt( reg2 ) * ( cos( b[7] ) + I * sin( b[7] ) ) ; 
    // cheaper method of reconstruction used, I have performed the matrix
    // multiplication by hand and simplified a little
    a[0] = a1 ; 
    a[1] = b[0] + I * b[1] ; 
    a[2] = b[2] + I * b[3] ;
    register const double complex star_a0 = conj( a[ 0 ] ) ; 
    const double complex ca1 = c * a[1] ;
    const double complex ca2 = c * a[2] ;
    a[3] = ( b[4] + I * b[5] ) ; 
    a[4] = (-conj( ca2 ) - a[3] * star_a0 * a[1] ) * NORM ;
    a[5] = ( conj( ca1 ) - star_a0 * a[2] * a[3] ) * NORM ;
    a[6] = c ;
    a[7] = ( conj( a[3] * a[2] ) - star_a0 * ca1 ) * NORM ;
    a[8] = ( -conj( a[3] * a[1] ) - star_a0 * ca2 ) * NORM ;
  }
#elif NC == 2
  const double R = sqrt( 1.0 - ( b[1] * b[1] + b[2] * b[2] ) ) ;
  a[0] = R * ( cos( b[0] ) + I * sin( b[0] ) ) ; 
  a[1] = b[1] + I * b[2] ;
  a[2] = -conj( a[1] ) ;
  a[3] = conj( a[0] ) ;
#endif
  // apparently sometimes need a reunitarisation because my constraints on how accurately
  // the code keeps stuff SU(3) have tightened and we are losing some backwards
  // compatibility here, so I reunit on a read of each 8 parameter configuration
  return ;
}
#endif

// same as above but for the antihermitian matrices
void
rebuild_antihermitian( GLU_complex a[ NCNC ] , 
		       const GLU_complex b[ HERMSIZE ] )
{
#if NC == 3
  a[0] = b[0] ;
  a[1] = b[1] ;
  a[2] = b[2] ;
  a[3] = -conj( a[1] ) ;
  a[4] = b[3] ;
  a[5] = b[4] ;
  a[6] = -conj( a[2] ) ;
  a[7] = -conj( a[5] ) ;
  a[8] = -a[0]-a[4] ;
#elif NC == 2 
  a[0] = b[0] ;
  a[1] = b[1] ;
  a[2] = -conj( a[1] ) ;
  a[3] = -a[0] ;
#else
  GLU_complex res = 0.0 ;
  size_t i , j , k , idx = 0 ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = i ; j < NC ; j++ ) {
      // match up the conjugates
      a[ j + NC * i ] = b[ idx ] ;
      a[ i + NC * j ] = -conj( b[ idx ] ) ;
      idx ++ ;
    }
  }
  // OK so we complete the last term as minus the sum of the diagonal
  // assuring tracelessness
  for( k = 0 ; k < NC-1 ; k++ ) {
    res += a[ k*(NC+1) ] ;
  }
  a[ NCNC - 1 ] = -res ;  
#endif
  return ;
}

// rebuild the full hermitian matrix from the packed hermitian
void 
rebuild_hermitian( GLU_complex a[ NCNC ] , 
		   const GLU_complex b[ HERMSIZE ] )
{
#if NC == 3
  a[0] = b[0] ;
  a[1] = b[1] ;
  a[2] = b[2] ;
  a[3] = conj( a[1] ) ;
  a[4] = b[3] ;
  a[5] = b[4] ;
  a[6] = conj( a[2] ) ;
  a[7] = conj( a[5] ) ;
  a[8] = -a[0]-a[4] ;
#elif NC == 2 
  a[0] = b[0] ;
  a[1] = b[1] ;
  a[2] = conj( a[1] ) ;
  a[3] = -a[0] ;
#else
  GLU_complex res = 0.0 ;
  int i , j , k , idx = 0 ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = i ; j < NC ; j++ ) {
      // match up the conjugates
      a[ j + NC * i ] = b[ idx ] ;
      a[ i + NC * j ] = conj( b[ idx ] ) ;
      idx ++ ;
    }
  }
  // OK so we complete the last term as minus the sum of the diagonal
  // assuring tracelessness
  for( k = 0 ; k < NC-1 ; k++ ) {
    res += a[ k*(NC+1) ] ;
  }
  a[ NCNC - 1 ] = -res ;  
#endif
  return ;
}

// shortening of the gauge matrices
#if (NC<4)
void
shorten( GLU_real *__restrict a ,
	 const GLU_complex b[ NCNC ] )
{
#if NC==3
  a[0] = creal( b[1] ) ; a[1] = cimag( b[1] ) ; 
  a[2] = creal( b[2] ) ; a[3] = cimag( b[2] ) ; 
  a[4] = creal( b[3] ) ; a[5] = cimag( b[3] ) ; 
  a[6] = carg( b[0] )  ; 
  a[7] = carg( b[6] )  ; 
#elif NC==2
  a[0] = carg( b[0] ) ;
  a[1] = creal( b[1] ) ;
  a[2] = cimag( b[1] ) ;
#endif
  return ;
}
#endif

// Determinant passed by reference
void
speed_det( GLU_complex *__restrict dt ,
	   const GLU_complex U[ NCNC ] )
{
#if NC == 3
  *dt = U[0] * ( U[4] * U[8] - U[5] * U[7] ) -	\
    U[1] * ( U[3] * U[8] - U[5] * U[6] ) +	\
    U[2] * ( U[3] * U[7] - U[4] * U[6] ) ; 
#elif NC == 2
  *dt = U[0] * U[3] - U[1] * U[2] ; 
#else
  *dt = (GLU_complex)LU_det( NC , U ) ;
#endif
  return ;
}

// Trace 
void
speed_trace( GLU_complex *res ,
	     const GLU_complex U[ NCNC ] ) 
{
#if NC == 3
  *res = U[0] + U[4] + U[8] ; 
#elif NC == 2
  *res = ( U[0] + U[3] )  ;
#else
  const GLU_real *pU = (const GLU_real*)U ;
  register double sumr = 0.0 , sumi = 0.0 ;
  size_t i ;
  for( i = 0 ; i < NC ; i++ ) { 
    sumr += ( *(pU++) ) ;
    sumi += ( *(pU++) ) ;
    pU += 2*NC ;
  }
  *res = sumr + I * sumi ;
#endif
  return ;
}

// real part of the trace by reference in *res
void
speed_trace_Re( double *__restrict res ,
		const GLU_complex U[ NCNC ] ) 
{
#if NC == 3
  *res = (double)creal( U[0] ) + (double)creal( U[4] ) + 
         (double)creal( U[8] ) ; 
#elif NC == 2
  *res = (double)creal( U[0] ) + (double)creal( U[3] ) ;
#else
  const GLU_real *pU = (const GLU_real*)U ;
  register double sumr = 0.0 ;
  size_t i ;
  for( i = 0 ; i < NC ; i++ ) { 
    sumr += (double)( *(pU) ) ;
    pU += 2*( NC + 1 ) ;
  }
  *res = sumr ;
#endif
  return ;
}

// Trace //
GLU_complex
trace( const GLU_complex U[ NCNC ] )
{
#if NC == 3
  return U[0] + U[4] + U[8] ;  
#elif NC == 2
  return U[0] + U[3] ;
#else
  const GLU_real *pU = (const GLU_real*)U ;
  register double sumr = 0.0 , sumi = 0.0 ;
  size_t i ;
  for( i = 0 ; i < NC ; i++ ) { 
    sumr += ( *(pU++) ) ;
    sumi += ( *(pU++) ) ;
    pU += 2*NC ; 
  }
  return sumr + I * sumi ;
#endif
}

// Trace of the product of two matrices //
void
trace_ab( GLU_complex *__restrict tr , 
	  const GLU_complex a[ NCNC ] , 
	  const GLU_complex b[ NCNC ] )
{
#if NC == 3
  *tr = a[0] * b[0] + a[1] * b[3] + a[2] * b[6] +	\
    a[3] * b[1] + a[4] * b[4] + a[5] * b[7] +		\
    a[6] * b[2] + a[7] * b[5] + a[8] * b[8] ;
#elif NC == 2
  *tr = a[0] * b[0] + a[2] * b[1] +\
    a[1] * b[2] + a[3] * b[3] ;
#else
  register GLU_real sumr = 0.0 , sumi = 0.0 ;
  size_t i, j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sumr += creal( a[ j + NC * i ] ) * creal( b[ i + NC * j ] ) -
	cimag( a[ j + NC * i ] ) * cimag( b[ i + NC * j ] ) ;
      sumi += creal( a[ j + NC * i ] ) * cimag( b[ i + NC * j ] ) +
	cimag( a[ j + NC * i ] ) * creal( b[ i + NC * j ] ) ;
    }
  }
  *tr = sumr + I * sumi ;
#endif
  return ;
}

// Trace of the product of two matrices ( b being daggered ) //
void
trace_ab_dag( GLU_complex *__restrict tr , 
	      const GLU_complex a[ NCNC ] , 
	      const GLU_complex b[ NCNC ] )
{
#if NC == 3
  *tr = a[0] * conj( b[0] ) + a[1] * conj( b[1] ) + a[2] * conj( b[2] ) + \
    a[3] * conj( b[3] ) + a[4] * conj( b[4] ) + a[5] * conj( b[5] ) +	\
    a[6] * conj( b[6] ) + a[7] * conj( b[7] ) + a[8] * conj( b[8] ) ;    
#elif NC == 2
  *tr = a[0] * conj( b[0] ) + a[1] * conj( b[1] ) + a[2] * conj( b[2] ) + \
    a[3] * conj( b[3] ) ;
#else
  register GLU_real sumr = 0.0 , sumi = 0.0 ;
  size_t i ;
  #if NC%2==0
  for( i = 0 ; i < NCNC ; i+=2 ) {
    sumr += creal( a[i] ) * creal( b[i] ) + cimag( a[i] ) * cimag( b[i] ) ;
    sumi += creal( a[i] ) * cimag( b[i] ) - creal( b[i] ) * cimag( a[i] ) ;
    sumr += creal( a[i+1] ) * creal( b[i+1] ) + cimag( a[i+1] ) * cimag( b[i+1] ) ;
    sumi += creal( a[i+1] ) * cimag( b[i+1] ) - creal( b[i+1] ) * cimag( a[i+1] ) ;
  }
  #else
  for( i = 0 ; i < NCNC ; i++ ) {
    sumr += creal( a[i] ) * creal( b[i] ) + cimag( a[i] ) * cimag( b[i] ) ;
    sumi += creal( a[i] ) * cimag( b[i] ) - creal( b[i] ) * cimag( a[i] ) ;
  }
  #endif
  *tr = sumr - I * sumi ;
#endif
  return ;
}

// Trace of the product of two matrices ( b being daggered ) //
void
trace_ab_dag_Re( GLU_real*__restrict tr , 
		 const GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] )
{
#if NC == 3
  *tr  = creal( a[0] ) * creal( b[0] ) + cimag( a[0] ) * cimag( b[0] ) ;
  *tr += creal( a[1] ) * creal( b[1] ) + cimag( a[1] ) * cimag( b[1] ) ;
  *tr += creal( a[2] ) * creal( b[2] ) + cimag( a[2] ) * cimag( b[2] ) ;
  *tr += creal( a[3] ) * creal( b[3] ) + cimag( a[3] ) * cimag( b[3] ) ;
  *tr += creal( a[4] ) * creal( b[4] ) + cimag( a[4] ) * cimag( b[4] ) ;
  *tr += creal( a[5] ) * creal( b[5] ) + cimag( a[5] ) * cimag( b[5] ) ;
  *tr += creal( a[6] ) * creal( b[6] ) + cimag( a[6] ) * cimag( b[6] ) ;
  *tr += creal( a[7] ) * creal( b[7] ) + cimag( a[7] ) * cimag( b[7] ) ;
  *tr += creal( a[8] ) * creal( b[8] ) + cimag( a[8] ) * cimag( b[8] ) ;
#elif NC == 2
  *tr = creal( a[0] ) * creal( b[0] ) + cimag( a[0] ) * cimag( b[0] ) ;
  *tr += creal( a[1] ) * creal( b[1] ) + cimag( a[1] ) * cimag( b[1] ) ;
  *tr += creal( a[2] ) * creal( b[2] ) + cimag( a[2] ) * cimag( b[2] ) ;
  *tr += creal( a[3] ) * creal( b[3] ) + cimag( a[3] ) * cimag( b[3] ) ;
#else
  register GLU_real sumr = 0.0 ;
  size_t i ;
  #if NC%2==0
  for( i = 0 ; i < NCNC ; i+=2 ) {
    sumr += creal( a[i] ) * creal( b[i] ) + cimag( a[i] ) * cimag( b[i] ) ;
    sumr += creal( a[i+1] ) * creal( b[i+1] ) + cimag( a[i+1] ) * cimag( b[i+1] ) ;
  }
  #else
  for( i = 0 ; i < NCNC ; i++ ) {
    sumr += creal( a[i] ) * creal( b[i] ) + cimag( a[i] ) * cimag( b[i] ) ;
  }
  #endif
  *tr = sumr ;
#endif
  return ;
}

// Trace of the product of two (different!) hermitian matrices is always real //
void
trace_ab_herm( GLU_real *__restrict tr , 
	       const GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC ] )
{
#if NC == 3
  // this is an identity
  *tr = creal( a[0] ) * creal( b[0] ) + creal( a[4] ) * creal( b[4] ) + creal( a[8] ) * creal( b[8] ) +
    2.0 * ( creal( a[1] ) * creal( b[1] ) + cimag( a[1] ) * cimag( b[1] ) + \
	    creal( a[2] ) * creal( b[2] ) + cimag( a[2] ) * cimag( b[2] ) + \
	    creal( a[5] ) * creal( b[5] ) + cimag( a[5] ) * cimag( b[5] ) ) ;
#elif NC == 2
  *tr = 2.0 * ( creal( a[0] ) * creal( b[0] ) +
		creal( a[1] ) * creal( b[1] ) + cimag( a[1] ) * cimag( b[1] ) ) ;
#else
  // loop over upper triangular taking the product
  *tr = 0.0 ;
  register double sum = 0.0 ;
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    sum += creal( a[ i + NC *i ] ) * creal( b[ i + NC *i ] ) ;
    for( j = i + 1 ; j < NC ; j++ ) {
      sum += 2.0 * ( creal( a[ j + NC *i ] ) * creal( b[ j + NC *i ] ) +
		     cimag( a[ j + NC *i ] ) * cimag( b[ j + NC *i ] ) ) ;
    }
  }
  *tr = (GLU_real)sum ;
#endif
  return ;
}

// Trace of the product of two shortened Hermitian matrices //
void
trace_ab_herm_short( GLU_real *__restrict tr , 
		     const GLU_complex a[ HERMSIZE ] , 
		     const GLU_complex b[ HERMSIZE ] )
{
#if NC == 3
  *tr = ( creal( a[0] ) * creal( b[3] ) + creal( a[3] ) * creal( b[0] ) ) + \
    2.0 * ( creal( a[0] ) * creal( b[0] ) + creal( a[3] ) * creal( b[3] ) + \
	    creal( a[1] ) * creal( b[1] ) + cimag( a[1] ) * cimag( b[1] ) + \
	    creal( a[2] ) * creal( b[2] ) + cimag( a[2] ) * cimag( b[2] ) + \
	    creal( a[4] ) * creal( b[4] ) + cimag( a[4] ) * cimag( b[4] ) ) ;
#elif NC == 2
  *tr = 2.0 * ( creal( a[0] ) * creal( b[0] ) +
		creal( a[1] ) * creal( b[1] ) + cimag( a[1] ) * cimag( b[1] ) ) ;
#else
  // loop over upper triangular taking the product
  *tr = 0.0 ;
  register double suma = 0.0 , sumb = 0.0 , sum = 0.0 ;
  size_t i , j , idx = 0 ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    sum += creal( a[idx] ) * creal( b[idx] ) ;
    suma += creal( a[idx] ) ;
    sumb += creal( b[idx] ) ;
    idx++ ;
    for( j = i + 1 ; j < NC ; j++ ) {
      sum += 2.0 * ( creal( a[idx] ) * creal( b[idx] )
		     + cimag( a[idx] ) * cimag( b[idx] ) ) ;
      idx++ ;
    }
  }
  *tr = (GLU_real)( sum + suma * sumb ) ;
#endif
}

//Trace of the product of the square of two ( hermitian ) matrices //
void
trace_prod_herm( GLU_real *__restrict tr ,
		 const GLU_complex a[ NCNC ] )
{  
#if NC == 3
  *tr = creal( a[0] ) * creal( a[0] )					\
    + creal( a[4] ) * creal( a[4] )					\
    + creal( a[8] ) * creal( a[8] )					\
    + 2.0 * ( creal( a[1] ) * creal( a[1] ) + cimag( a[1] ) * cimag( a[1] ) \
	      + creal( a[2] ) * creal( a[2] ) + cimag( a[2] ) * cimag( a[2] ) \
	      + creal( a[5] ) * creal( a[5] ) + cimag( a[5] ) * cimag( a[5] ) ) ;
#elif NC == 2
  *tr = 2.0 * ( creal( a[0] ) * creal( a[0] ) +				\
		creal( a[1] ) * creal( a[1] ) + cimag( a[1] ) * cimag( a[1] ) ) ;
#else
  // remove the modulus, and pointer-ise gave me 2x maybe a little compiler
  // dependent
  const GLU_complex *pA = a ;
  register double sum = 0.0 ;
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    sum += creal( *pA ) * creal( *pA ) ;
    pA++ ;
    for( j = i + 1 ; j < NC ; j++ ) {
      sum += 2.0 * ( creal( *pA ) * creal( *pA ) 
		     + cimag( *pA ) * cimag( *pA ) ) ;
      pA++ ;
    }
    pA += i+1 ;
  }
  *tr = sum ;
#endif 
  return ;
}

// transpose some matrix b into a //
void
transpose( GLU_complex a[ NCNC ] , 
	   const GLU_complex b[ NCNC ] )
{
#if NC == 3
  a[0] = b[0] ; a[1] = b[3] ; a[2] = b[6] ; 
  a[3] = b[1] ; a[4] = b[4] ; a[5] = b[7] ; 
  a[6] = b[2] ; a[7] = b[5] ; a[8] = b[8] ; 
#elif NC == 2
  a[0] = b[0] ; a[1] = b[2] ; 
  a[2] = b[1] ; a[3] = b[3] ;
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      a[ j + i * NC ] = b[ i + j * NC ] ;
    }
  }
#endif
  return ;
}

// This function prints to the terminal the visualization of our matrices //
void 
write_matrix( const GLU_complex U[ NCNC ] )
{
  size_t i ; 
  for( i = 0 ; i < NCNC ; i++ ) {
    if( i % NC == 0 ) {
      fprintf( stdout , "\n" ) ;
    } 
    fprintf( stdout , "%f %f " , creal( U[i] ) , cimag( U[i] ) ) ; 
  }
  fprintf( stdout , "\n \n" ) ; 
  return ;
}

// This function prints to the terminal the array in cform so it can be checked in separate codes//
void 
write_matrix_cform( const GLU_complex U[ NCNC ] )
{
  size_t i ; 
  fprintf( stdout , "\n{") ;
  for( i = 0 ; i < NCNC ; i++ ) {
    fprintf( stdout , " %f + I*%f " , creal( U[i] ) , cimag( U[i] ) ) ;
    if( i < NCNC-1 ){ fprintf( stdout , "," ) ; }
  }
  fprintf( stdout , "} ;\n" ) ;
  return ;
}

// This function prints to the terminal the Mathematica form visualization of our matrices //
void 
write_matrix_mathematica( const GLU_complex U[ NCNC ] ) 
{
  size_t i ; 
  fprintf( stdout , "\n{{" ) ; 
  for( i = 0 ; i < NCNC ; i++ ) {
    if( ( i%NC == 0 ) && ( i != 0 ) ) {
      fprintf( stdout , "} , " ) ;
      fprintf( stdout , "{" ) ; 
    } else {
      if( i != 0 ) {
	fprintf( stdout , " , " ) ;
      }
    }
    fprintf( stdout , "%1.15f+I*%1.15f" , creal( U[i] ) , cimag( U[i] ) ) ; 
  }
  fprintf( stdout , "}}\n \n" ) ; 
  return ;
}

// zero a matrix //
void
zero_mat( GLU_complex a[ NCNC ] ) 
{
  // should probably just be a call to memset
#if NC == 3
  a[0] = 0. ;  a[1] = 0. ;  a[2] = 0. ; 
  a[3] = 0. ;  a[4] = 0. ;  a[5] = 0. ; 
  a[6] = 0. ;  a[7] = 0. ;  a[8] = 0. ; 
#elif NC == 2
  a[0] = 0. ;  a[1] = 0. ;  
  a[2] = 0. ;  a[3] = 0. ;  
#else
  // should call to memset here?
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = 0. ;
  }
#endif
  return ;
}





