/*
    Copyright 2013 Renwick James Hudspith

    This file (MMUL.c) is part of GLU.

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
   @file MMUL.c
   @brief matrix-matrix multiply
 */

#include "Mainfile.h"

#if !( defined HAVE_IMMINTRIN_H ) || ( defined SINGLE_PREC )

#ifndef multab
// simple matrix multiplication
void 
multab( GLU_complex a[ NCNC ] , 
	const GLU_complex b[ NCNC ] , 
	const GLU_complex c[ NCNC ] )
{
#if NC==3
  a[0] = b[0] * c[0] + b[1] * c[3] + b[2] * c[6] ;	\
  a[1] = b[0] * c[1] + b[1] * c[4] + b[2] * c[7] ;	\
  a[2] = b[0] * c[2] + b[1] * c[5] + b[2] * c[8] ;	\
  a[3] = b[3] * c[0] + b[4] * c[3] + b[5] * c[6] ;	\
  a[4] = b[3] * c[1] + b[4] * c[4] + b[5] * c[7] ;	\
  a[5] = b[3] * c[2] + b[4] * c[5] + b[5] * c[8] ;	\
  a[6] = b[6] * c[0] + b[7] * c[3] + b[8] * c[6] ;	\
  a[7] = b[6] * c[1] + b[7] * c[4] + b[8] * c[7] ;	\
  a[8] = b[6] * c[2] + b[7] * c[5] + b[8] * c[8] ;	
#elif NC==2
  a[0] = b[0] * c[0] + b[1] * c[2] ;		\
  a[1] = b[0] * c[1] + b[1] * c[3] ;		\
  a[2] = b[2] * c[0] + b[3] * c[2] ;		\
  a[3] = b[2] * c[1] + b[3] * c[3] ;		
#else
  // slow and stupid version
  size_t i , j , m ;
  register GLU_complex sum ;
  register GLU_real REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = 0.0 ;
      for( m = 0 ; m < NC ; m++  ) {
	REB = creal( b[ m + NC*i ] ) ; IMB = cimag( b[ m + NC*i ] ) ;
	REC = creal( c[ j + m*NC ] ) ; IMC = cimag( c[ j + m*NC ] ) ;
	sum += REB * REC - IMB * IMC + I * ( REB * IMC + IMB * REC ) ; 
      }
      a[ j + NC*i ] = sum ;
    }
  }
#endif
  return ;
}
#endif

#ifndef multab_atomic_left
// simple matrix multiplication ( left multiply ) a = b * a 
void 
multab_atomic_left( GLU_complex a[ NCNC ] , 
		    const GLU_complex b[ NCNC ] )
{
#if NC == 3
  GLU_complex C0 = b[0] * a[0] + b[1] * a[3] + b[2] * a[6] ;	\
  GLU_complex C1 = b[3] * a[0] + b[4] * a[3] + b[5] * a[6] ;	\
  GLU_complex C2 = b[6] * a[0] + b[7] * a[3] + b[8] * a[6] ;	\
  a[0] = C0 ; a[3] = C1 ; a[6] = C2 ;				\
  C0 = b[0] * a[1] + b[1] * a[4] + b[2] * a[7] ;		\
  C1 = b[3] * a[1] + b[4] * a[4] + b[5] * a[7] ;		\
  C2 = b[6] * a[1] + b[7] * a[4] + b[8] * a[7] ;		\
  a[1] = C0 ; a[4] = C1 ; a[7] = C2 ;				\
  C0 = b[0] * a[2] + b[1] * a[5] + b[2] * a[8] ;		\
  C1 = b[3] * a[2] + b[4] * a[5] + b[5] * a[8] ;		\
  C2 = b[6] * a[2] + b[7] * a[5] + b[8] * a[8] ;		\
  a[2] = C0 ; a[5] = C1 ; a[8] = C2 ;
#elif NC == 2
  GLU_complex C0 = b[0] * a[0] + b[1] * a[2] ;			\
  GLU_complex C1 = b[2] * a[0] + b[3] * a[2] ;			\
  a[0] = C0 ; a[2] = C1 ;					\
  C0 = b[0] * a[1] + b[1] * a[3] ;				\
  C1 = b[2] * a[1] + b[3] * a[3] ;				\
  a[1] = C0 ; a[3] = C1 ;					
#else
  // standard looped version
  size_t i , j , m ;
  GLU_complex R[ NC ] ;
  register GLU_complex sum ;
  for( i = 0 ; i < NC ; i++ ) { // loop cols
    for( j = 0 ; j < NC ; j ++ ) { // loop rows
      sum = 0.0 ;
      for( m = 0 ; m < NC ; m ++  ) { // loop elements in row or column
        sum += b[ m + j*NC ] * a[ i + m*NC ] ;
      }	
      R[j] = sum ;
    }
    // copy back over to a ...
    for( m = 0 ; m < NC ; m++ ) {
      a[ i + m*NC ] = R[m] ;
    }
  }
#endif
  return ;
}
#endif

#ifndef multab_atomic
// simple matrix multiplication a = a * b
void 
multab_atomic_right( GLU_complex a[ NCNC ] , 
		     const GLU_complex b[ NCNC ] )
{
#if NC==3
  GLU_complex R0 = a[0] * b[0] + a[1] * b[3] + a[2] * b[6] ;	\
  GLU_complex R1 = a[0] * b[1] + a[1] * b[4] + a[2] * b[7] ;	\
  GLU_complex R2 = a[0] * b[2] + a[1] * b[5] + a[2] * b[8] ;	\
  a[0] = R0 ; a[1] = R1 ; a[2] = R2 ;				\
  R0 = a[3] * b[0] + a[4] * b[3] + a[5] * b[6] ;		\
  R1 = a[3] * b[1] + a[4] * b[4] + a[5] * b[7] ;		\
  R2 = a[3] * b[2] + a[4] * b[5] + a[5] * b[8] ;		\
  a[3] = R0 ; a[4] = R1 ; a[5] = R2 ;				\
  R0 = a[6] * b[0] + a[7] * b[3] + a[8] * b[6] ;		\
  R1 = a[6] * b[1] + a[7] * b[4] + a[8] * b[7] ;		\
  R2 = a[6] * b[2] + a[7] * b[5] + a[8] * b[8] ;		\
  a[6] = R0 ; a[7] = R1 ; a[8] = R2 ;				
#elif NC==2
  GLU_complex R0 = a[0] * b[0] + a[1] * b[2] ;		\
  GLU_complex R1 = a[0] * b[1] + a[1] * b[3] ;		\
  a[0] = R0 ; a[1] = R1 ;				\
  R0 = a[2] * b[0] + a[3] * b[2] ;			\
  R1 = a[2] * b[1] + a[3] * b[3] ;			\
  a[2] = R0 ; a[3] = R1 ;				
#else
  // slow and stupid loopy version
  size_t i , j , m ;
  GLU_complex R[ NC ] ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j ++ ) {
      R[j] = 0.0 ;
      for( m = 0 ; m < NC ; m ++  ) {
        R[ j ] += a[ m + i*NC ] * b[ j + m*NC ] ;
      }	
    }
    // copy back over to a ...
    for( j = 0 ; j < NC ; j ++ ) {
      a[ j + i*NC ] = R[j] ;
    }
  }
#endif
  return ;
}
#endif

#ifndef multab_suNC

void //__attribute__((hot))
multab_suNC( GLU_complex a[ NCNC ] , 
	     const GLU_complex b[ NCNC ] , 
	     const GLU_complex c[ NCNC ] )
{
#if NC == 3
  //a[0] = b[0] * c[0] + b[1] * c[3] + b[2] * c[6] ; 
  a[0] = creal( b[0] ) * creal( c[0] ) - cimag( b[0] ) * cimag( c[0] ) + I * ( creal( b[0] ) * cimag( c[0] ) + creal( c[0] ) * cimag( b[0] ) ) \
       + creal( b[1] ) * creal( c[3] ) - cimag( b[1] ) * cimag( c[3] ) + I * ( creal( b[1] ) * cimag( c[3] ) + creal( c[3] ) * cimag( b[1] ) ) \
       + creal( b[2] ) * creal( c[6] ) - cimag( b[2] ) * cimag( c[6] ) + I * ( creal( b[2] ) * cimag( c[6] ) + creal( c[6] ) * cimag( b[2] ) ) ; 
  //a[1] = b[0] * c[1] + b[1] * c[4] + b[2] * c[7] ; 
  a[1] = creal( b[0] ) * creal( c[1] ) - cimag( b[0] ) * cimag( c[1] ) + I * ( creal( b[0] ) * cimag( c[1] ) + creal( c[1] ) * cimag( b[0] ) ) \
       + creal( b[1] ) * creal( c[4] ) - cimag( b[1] ) * cimag( c[4] ) + I * ( creal( b[1] ) * cimag( c[4] ) + creal( c[4] ) * cimag( b[1] ) ) \
       + creal( b[2] ) * creal( c[7] ) - cimag( b[2] ) * cimag( c[7] ) + I * ( creal( b[2] ) * cimag( c[7] ) + creal( c[7] ) * cimag( b[2] ) ) ; 
  //a[2] = b[0] * c[2] + b[1] * c[5] + b[2] * c[8] ; 
  a[2] = creal( b[0] ) * creal( c[2] ) - cimag( b[0] ) * cimag( c[2] ) + I * ( creal( b[0] ) * cimag( c[2] ) + creal( c[2] ) * cimag( b[0] ) ) \
       + creal( b[1] ) * creal( c[5] ) - cimag( b[1] ) * cimag( c[5] ) + I * ( creal( b[1] ) * cimag( c[5] ) + creal( c[5] ) * cimag( b[1] ) ) \
       + creal( b[2] ) * creal( c[8] ) - cimag( b[2] ) * cimag( c[8] ) + I * ( creal( b[2] ) * cimag( c[8] ) + creal( c[8] ) * cimag( b[2] ) ) ; 
  // middle row //
  //a[3] = b[3] * c[0] + b[4] * c[3] + b[5] * c[6] ; 
  a[3] = creal( b[3] ) * creal( c[0] ) - cimag( b[3] ) * cimag( c[0] ) + I * ( creal( b[3] ) * cimag( c[0] ) + creal( c[0] ) * cimag( b[3] ) ) \
       + creal( b[4] ) * creal( c[3] ) - cimag( b[4] ) * cimag( c[3] ) + I * ( creal( b[4] ) * cimag( c[3] ) + creal( c[3] ) * cimag( b[4] ) ) \
       + creal( b[5] ) * creal( c[6] ) - cimag( b[5] ) * cimag( c[6] ) + I * ( creal( b[5] ) * cimag( c[6] ) + creal( c[6] ) * cimag( b[5] ) ) ; 
  // a[4] = b[3] * c[1] + b[4] * c[4] + b[5] * c[7] ; 
  a[4] = creal( b[3] ) * creal( c[1] ) - cimag( b[3] ) * cimag( c[1] ) + I * ( creal( b[3] ) * cimag( c[1] ) + creal( c[1] ) * cimag( b[3] ) ) \
       + creal( b[4] ) * creal( c[4] ) - cimag( b[4] ) * cimag( c[4] ) + I * ( creal( b[4] ) * cimag( c[4] ) + creal( c[4] ) * cimag( b[4] ) ) \
       + creal( b[5] ) * creal( c[7] ) - cimag( b[5] ) * cimag( c[7] ) + I * ( creal( b[5] ) * cimag( c[7] ) + creal( c[7] ) * cimag( b[5] ) ) ; 
  //a[5] = b[3] * c[2] + b[4] * c[5] + b[5] * c[8] ; 
  a[5] = creal( b[3] ) * creal( c[2] ) - cimag( b[3] ) * cimag( c[2] ) + I * ( creal( b[3] ) * cimag( c[2] ) + creal( c[2] ) * cimag( b[3] ) ) \
       + creal( b[4] ) * creal( c[5] ) - cimag( b[4] ) * cimag( c[5] ) + I * ( creal( b[4] ) * cimag( c[5] ) + creal( c[5] ) * cimag( b[4] ) ) \
       + creal( b[5] ) * creal( c[8] ) - cimag( b[5] ) * cimag( c[8] ) + I * ( creal( b[5] ) * cimag( c[8] ) + creal( c[8] ) * cimag( b[5] ) ) ; 
  // bottom row // as a completion of the top two
  //a[6] = conj( a[1]  *  a[5] - a[2]  *  a[4] ) ; 
  a[6] = creal( a[1] ) * creal( a[5] ) - cimag( a[1] ) * cimag( a[5] ) - I * ( creal( a[1] ) * cimag( a[5] ) + creal( a[5] ) * cimag( a[1] ) ) \
       - creal( a[2] ) * creal( a[4] ) + cimag( a[2] ) * cimag( a[4] ) + I * ( creal( a[2] ) * cimag( a[4] ) + creal( a[4] ) * cimag( a[2] ) ) ; 
  // a[7] = conj( a[2]  *  a[3] - a[0]  *  a[5] ) ; 
  a[7] = creal( a[2] ) * creal( a[3] ) - cimag( a[2] ) * cimag( a[3] ) - I * ( creal( a[2] ) * cimag( a[3] ) + creal( a[3] ) * cimag( a[2] ) ) \
       - creal( a[0] ) * creal( a[5] ) + cimag( a[0] ) * cimag( a[5] ) + I * ( creal( a[0] ) * cimag( a[5] ) + creal( a[5] ) * cimag( a[0] ) ) ; 
  // a[8] = conj( a[0]  *  a[4] - a[1]  *  a[3] ) ; 
  a[8] = creal( a[0] ) * creal( a[4] ) - cimag( a[0] ) * cimag( a[4] ) - I * ( creal( a[0] ) * cimag( a[4] ) + creal( a[4] ) * cimag( a[0] ) ) \
       - creal( a[1] ) * creal( a[3] ) + cimag( a[1] ) * cimag( a[3] ) + I * ( creal( a[1] ) * cimag( a[3] ) + creal( a[3] ) * cimag( a[1] ) ) ;
#elif NC == 2
  a[0] = creal( b[0] ) * creal( c[0] ) - cimag( b[0] ) * cimag( c[0] ) + \
    creal( b[1] ) * creal( c[2] ) - cimag( b[1] ) * cimag( c[2] ) +	\
    I * ( creal( b[0] ) * cimag( c[0] ) + creal( c[0] ) * cimag( b[0] ) + \
	  creal( b[1] ) * cimag( c[2] ) + creal( c[2] ) * cimag( b[1] ) ) ;  

  a[1] = creal( b[0] ) * creal( c[1] ) - cimag( b[0] ) * cimag( c[1] ) + \
    creal( b[1] ) * creal( c[3] ) - cimag( b[1] ) * cimag( c[3] ) +	\
    I * ( creal( b[0] ) * cimag( c[1] ) + creal( c[1] ) * cimag( b[0] ) + \
	  creal( b[1] ) * cimag( c[3] ) + creal( c[3] ) * cimag( b[1] ) ) ; 

  //a[2] = b[2]*c[0] + b[3]*c[2]
  a[2] = -conj( a[1] ) ;  
   // a[3] = b[2]*c[1] + b[3]*c[3]
  a[3] = conj( a[0] ) ; 

#else
  return multab_suNC( a , b , c ) ;
#endif
  return ;
}

#endif

#endif // <immintrin.h>
