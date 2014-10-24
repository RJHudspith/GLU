/*
    Copyright 2013 Renwick James Hudspith

    This file (gramschmidt.c) is part of GLU.

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
   @file gramschmidt.c
   @brief reunitarisation procedures
 */

#include "Mainfile.h"

#include "GLU_rng.h"  // generate_NCxNC() is called

#if NC > 3

// gramschmidt projection V = V - V.U^{\dagger}
static void
project( V , U )
     GLU_complex *__restrict V ;
     const GLU_complex *__restrict U ;
{
  int i ;
  register double projRE = 0.0 , projIM = 0.0 ;
  for( i = 0 ; i < NC ; i ++ ) {
    projRE += ( double )( creal( V[i] ) * creal( U[i] ) ) ;
    projRE += ( double )( cimag( V[i] ) * cimag( U[i] ) ) ;
    projIM += ( double )( creal( U[i] ) * cimag( V[i] ) ) ;
    projIM -= ( double )( cimag( U[i] ) * creal( V[i] ) ) ; 
  } 
  for( i = 0 ; i < NC ; i++ ) {
    V[i] -= projRE * creal( U[i] ) - projIM * cimag( U[i] ) +
      I * ( projRE * cimag( U[i] ) + projIM * creal( U[i] ) ) ;
  }
  return ;
}

// normalize a vector //
INLINE_STATIC_VOID
vect_norm2( a )
     GLU_complex *__restrict a ;
{
  int mu ;
  register double norm = 0.0 ;
  for( mu = 0 ; mu < NC ; mu++ ) {
    norm += (double)( creal( a[mu] ) * creal( a[mu] ) ) +	\
            (double)( cimag( a[mu] ) * cimag( a[mu] ) ) ;
  }
  norm = 1. / sqrt( norm ) ;
  for( mu = 0 ; mu < NC ; mu++ ) {
    a[ mu ] *= norm ;
  }
  return ;
}
#endif

//reunitarise an SU(3) matrix sped up for our requirements
void 
reunit2( GLU_complex *__restrict U)
{
#if NC==3

  GLU_real *uu = (GLU_real*)U ;
  double v00R = *( uu + 0 ) ;
  double v00I = *( uu + 1 ) ;
  double v10R = *( uu + 2 ) ;
  double v10I = *( uu + 3 ) ;
  double v01R = *( uu + 6 ) ;
  double v01I = *( uu + 7 ) ;
  double v11R = *( uu + 8 ) ;
  double v11I = *( uu + 9 ) ;
  double v02R = *( uu + 12 ) ;
  double v02I = *( uu + 13 ) ;
  double v12R = *( uu + 14 ) ;
  double v12I = *( uu + 15 ) ;

  double v20R , v20I , v21R , v21I , v22R , v22I ;
  double cR , cI , norm ;
  
  //the procedure i use is a gram schmidt orthogonalization with a cross product
  //the first column is normalized and invariant second orthogonal third cross product
  norm = ( v00R * v00R + v00I * v00I ) ;
  norm = ( v01R * v01R + v01I * v01I ) + norm ;
  norm = ( v02R * v02R + v02I * v02I ) + norm ;	
  norm = 1. / sqrt( norm ) ;

  v00R *= norm ;
  v00I *= norm ;
  v01R *= norm ;
  v01I *= norm ;
  v02R *= norm ;
  v02I *= norm ;

  //initial gramschmidt
  cR = ( v10R * v00R ) + ( v10I * v00I ) ;
  cR = ( v11R * v01R ) + ( v12R * v02R ) + cR ;
  cR = ( v12I * v02I ) + ( v11I * v01I ) + cR ;

  cI = ( v00R * v10I - v10R * v00I ) ;
  cI = ( v01R * v11I - v11R * v01I ) + cI ;
  cI = ( v02R * v12I - v12R * v02I ) + cI ;

  //now with 25% more loop unrolling
  v10R -= ( cR * v00R - cI * v00I ) ; 
  v10I -= ( cR * v00I + cI * v00R ) ; 
  v11R -= ( cR * v01R - cI * v01I ) ; 
  v11I -= ( cR * v01I + cI * v01R ) ; 
  v12R -= ( cR * v02R - cI * v02I ) ; 
  v12I -= ( cR * v02I + cI * v02R ) ; 

  norm = v10R * v10R + v10I * v10I ;
  norm = v11R * v11R + v11I * v11I + norm ;
  norm = v12R * v12R + v12I * v12I + norm ;	
  norm = 1. / sqrt( norm ) ;

  v10R *= norm ;
  v10I *= norm ;
  v11R *= norm ;
  v11I *= norm ;
  v12R *= norm ;
  v12I *= norm ;

  //complete
  v20R = v01R * v12R - v01I * v12I ; 
  v20R = v20R - v02R * v11R + v02I * v11I ;
  v20I = v01R * v12I + v12R * v01I ;
  v20I = v02R * v11I + v11R * v02I - v20I ;

  v21R = v02R * v10R - v02I * v10I ; 
  v21R = v21R - v12R * v00R + v12I * v00I ;
  v21I = v02R * v10I + v10R * v02I ;
  v21I = v00R * v12I + v12R * v00I - v21I ;

  v22R = v00R * v11R - v00I * v11I ; 
  v22R = v22R - v01R * v10R + v01I * v10I ;
  v22I = v00R * v11I + v11R * v00I ;
  v22I = v01R * v10I + v10R * v01I - v22I ;
 
  *( ( (GLU_real*) U + 0 ) ) = v00R ;
  *( ( (GLU_real*) U + 1 ) ) = v00I ;
  *( ( (GLU_real*) U + 2 ) ) = v10R ;
  *( ( (GLU_real*) U + 3 ) ) = v10I ;
  *( ( (GLU_real*) U + 4 ) ) = v20R ;
  *( ( (GLU_real*) U + 5 ) ) = v20I ;
  *( ( (GLU_real*) U + 6 ) ) = v01R ;
  *( ( (GLU_real*) U + 7 ) ) = v01I ;
  *( ( (GLU_real*) U + 8 ) ) = v11R ;
  *( ( (GLU_real*) U + 9 ) ) = v11I ;
  *( ( (GLU_real*) U + 10 ) ) = v21R ;
  *( ( (GLU_real*) U + 11 ) ) = v21I ;
  *( ( (GLU_real*) U + 12 ) ) = v02R ;
  *( ( (GLU_real*) U + 13 ) ) = v02I ;
  *( ( (GLU_real*) U + 14 ) ) = v12R ;
  *( ( (GLU_real*) U + 15 ) ) = v12I ;
  *( ( (GLU_real*) U + 16 ) ) = v22R ;
  *( ( (GLU_real*) U + 17 ) ) = v22I ;
 
#elif NC == 2

  GLU_real *uu = (GLU_real*)U ;
  const double v00R = *( uu + 0 ) ;
  const double v00I = *( uu + 1 ) ;
  const double v01R = *( uu + 4 ) ;
  const double v01I = *( uu + 5 ) ;

  double norm = sqrt( v00R * v00R + v00I * v00I + v01R * v01R + v01I * v01I ) ;
  norm = 1. / norm ;

  U[0] = ( v00R * norm + I * v00I * norm ); 
  U[1] = -( v01R * norm - I * v01I * norm ) ; 
  U[2] = ( v01R * norm + I * v01I * norm ) ; 
  U[3] = ( v00R * norm - I * v00I * norm ); 

#else
  // need a modified gram-schmidt process, leaves the
  // bottom row untouched
  int i , j ;
  GLU_complex v[ NC - 1 ][ NC ] ;

  // equate our vector to our matrix , parallelism here is a bit extreme
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      v[i][j] = U[ j + ( NC * i ) ] ;
    }
  }

  // perform the modified gram schmidt
  int k ;
  for( k = 0 ; k < NC-1 ; k++ ) {
    for( i = 0 ; i < k ; i++ ) { 
      project( v[k] , v[i] ) ;      
    }
    // normalize
    vect_norm2( v[k] ) ;
    // and put back into U, need the index i to make sure the vector gets normalised
    for( j = 0 ; j < NC ; j++ ) {
      U[ j + ( NC * i ) ] = v[i][j] ;
    }
  }

  // OK so this process leaves the bottom row as was, need to complete
  GLU_complex *array = malloc( ( NC - 1 ) * ( NC - 1 ) * sizeof( GLU_complex ) ) ;
  for( i = (NCNC-NC) ; i < NCNC ; i++ ) { // our bona-fide minor index loops the bottom row
    int idx = 0 ;
    for( j = 0 ; j < ( NCNC - NC ) ; j++ ) {
      if( ( j%NC != i%NC ) ) { // remove columns  ( j/NC != i/NC ) is implicit!!
	// pack array
	array[idx++] = U[j] ;
      } 
    }
    // compute the determinant
    #if ( NC%2 == 0 )
    register const GLU_real mulfact = ( i % 2 == 0 ) ? -1.0 : 1.0 ; 
    #else
    register const GLU_real mulfact = ( i % 2 == 0 ) ? 1.0 : -1.0 ; 
    #endif
    // compute the determinant of the minors using the LU decomp
    U[i] = conj( mulfact * (GLU_complex)LU_det( NC-1 , array ) ) ;
  }
  free( array ) ;
#endif
  return ;
}

//reunitarise an SU(Nc) matrix
void 
reunit( GLU_complex Z[ NCNC ] ,
	const GLU_complex U[ NCNC ] ) 
{
  equiv( Z , U ) ;
  reunit2( Z ) ;
  return ;
}

// generate a random SU(N) matrix
void 
Sunitary_gen( GLU_complex Z[ NCNC ] )
{
  generate_NCxNC( Z ) ; // generate gaussian distributed elements of matrix
  reunit2( Z ) ;
  while( !is_unitary( Z ) ) {
    printf("not unitary! Redoing ") ;
    generate_NCxNC( Z ) ; 
    reunit2( Z ) ;
  }
  return ;
}





