/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (MMULdagdag_SUNC.c) is part of GLU.

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
   @file MMULdagdag_SUNC.c
   @brief computes \f$ a = b^{\dagger} c^{\dagger} \f$ matrix multiplication for specifically SU(#NC) matrices
 */
#include "Mainfile.h"

#if !( defined HAVE_IMMINTRIN_H ) || ( defined SINGLE_PREC )

#ifndef multab_dagdag_suNC

void //__attribute__((hot))
multab_dagdag_suNC( GLU_complex a[ NCNC ] , 
		    const GLU_complex b[ NCNC ] , 
		    const GLU_complex c[ NCNC ] )
{
#if NC == 3
  // top row //
  //a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] ) + conj( b[6] ) * conj( c[2] ) ; 
  a[0] = creal( b[0] ) * creal( c[0] ) - cimag( b[0] ) * cimag( c[0] ) - I * ( creal( b[0] ) * cimag( c[0] ) + creal( c[0] ) * cimag( b[0] ) ) \
       + creal( b[3] ) * creal( c[1] ) - cimag( b[3] ) * cimag( c[1] ) - I * ( creal( b[3] ) * cimag( c[1] ) + creal( c[1] ) * cimag( b[3] ) ) \
       + creal( b[6] ) * creal( c[2] ) - cimag( b[6] ) * cimag( c[2] ) - I * ( creal( b[6] ) * cimag( c[2] ) + creal( c[2] ) * cimag( b[6] ) ) ; 

  //a[1] = conj( b[0] ) * conj( c[3] ) + conj( b[3] ) * conj( c[4] ) + conj( b[6] ) * conj( c[5] ) ; 
  a[1] = creal( b[0] ) * creal( c[3] ) - cimag( b[0] ) * cimag( c[3] ) - I * ( creal( b[0] ) * cimag( c[3] ) + creal( c[3] ) * cimag( b[0] ) ) \
       + creal( b[3] ) * creal( c[4] ) - cimag( b[3] ) * cimag( c[4] ) - I * ( creal( b[3] ) * cimag( c[4] ) + creal( c[4] ) * cimag( b[3] ) ) \
       + creal( b[6] ) * creal( c[5] ) - cimag( b[6] ) * cimag( c[5] ) - I * ( creal( b[6] ) * cimag( c[5] ) + creal( c[5] ) * cimag( b[6] ) ) ; 

  //a[2] = conj( b[0] ) * conj( c[6] ) + conj( b[3] ) * conj( c[7] ) + conj( b[6] ) * conj( c[8] ) ; 
  a[2] = creal( b[0] ) * creal( c[6] ) - cimag( b[0] ) * cimag( c[6] ) - I * ( creal( b[0] ) * cimag( c[6] ) + creal( c[6] ) * cimag( b[0] ) ) \
       + creal( b[3] ) * creal( c[7] ) - cimag( b[3] ) * cimag( c[7] ) - I * ( creal( b[3] ) * cimag( c[7] ) + creal( c[7] ) * cimag( b[3] ) ) \
       + creal( b[6] ) * creal( c[8] ) - cimag( b[6] ) * cimag( c[8] ) - I * ( creal( b[6] ) * cimag( c[8] ) + creal( c[8] ) * cimag( b[6] ) ) ; 

  // middle row //
  //a[3] = conj( b[1] ) * conj( c[0] ) + conj( b[4] ) * conj( c[1] ) + conj( b[7] ) * conj( c[2] ) ; 
  a[3] = creal( b[1] ) * creal( c[0] ) - cimag( b[1] ) * cimag( c[0] ) - I * ( creal( b[1] ) * cimag( c[0] ) + creal( c[0] ) * cimag( b[1] ) ) \
       + creal( b[4] ) * creal( c[1] ) - cimag( b[4] ) * cimag( c[1] ) - I * ( creal( b[4] ) * cimag( c[1] ) + creal( c[1] ) * cimag( b[4] ) ) \
       + creal( b[7] ) * creal( c[2] ) - cimag( b[7] ) * cimag( c[2] ) - I * ( creal( b[7] ) * cimag( c[2] ) + creal( c[2] ) * cimag( b[7] ) ) ; 

  //a[4] = conj( b[1] ) * conj( c[3] ) + conj( b[4] ) * conj( c[4] ) + conj( b[7] ) * conj( c[5] ) ; 
  a[4] = creal( b[1] ) * creal( c[3] ) - cimag( b[1] ) * cimag( c[3] ) - I * ( creal( b[1] ) * cimag( c[3] ) + creal( c[3] ) * cimag( b[1] ) ) \
       + creal( b[4] ) * creal( c[4] ) - cimag( b[4] ) * cimag( c[4] ) - I * ( creal( b[4] ) * cimag( c[4] ) + creal( c[4] ) * cimag( b[4] ) ) \
       + creal( b[7] ) * creal( c[5] ) - cimag( b[7] ) * cimag( c[5] ) - I * ( creal( b[7] ) * cimag( c[5] ) + creal( c[5] ) * cimag( b[7] ) ) ; 

  //a[5] = conj( b[1] ) * conj( c[6] ) + conj( b[4] ) * conj( c[7] ) + conj( b[7] ) * conj( c[8] ) ; 
  a[5] = creal( b[1] ) * creal( c[6] ) - cimag( b[1] ) * cimag( c[6] ) - I * ( creal( b[1] ) * cimag( c[6] ) + creal( c[6] ) * cimag( b[1] ) ) \
       + creal( b[4] ) * creal( c[7] ) - cimag( b[4] ) * cimag( c[7] ) - I * ( creal( b[4] ) * cimag( c[7] ) + creal( c[7] ) * cimag( b[4] ) ) \
       + creal( b[7] ) * creal( c[8] ) - cimag( b[7] ) * cimag( c[8] ) - I * ( creal( b[7] ) * cimag( c[8] ) + creal( c[8] ) * cimag( b[7] ) ) ; 

  //bottom row 
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
  //a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[2] ) * conj( c[1] )  ; 
  a[0] = creal( b[0] ) * creal( c[0] ) - cimag( b[0] ) * cimag( c[0] ) - I * ( creal( b[0] ) * cimag( c[0] ) + creal( c[0] ) * cimag( b[0] ) ) \
       + creal( b[2] ) * creal( c[1] ) - cimag( b[2] ) * cimag( c[1] ) - I * ( creal( b[2] ) * cimag( c[1] ) + creal( c[1] ) * cimag( b[2] ) ) ;  

  //a[1] = conj( b[0] ) * conj( c[2] ) + conj( b[2] ) * conj( c[3] )  ; 
  a[1] = creal( b[0] ) * creal( c[2] ) - cimag( b[0] ) * cimag( c[2] ) - I * ( creal( b[0] ) * cimag( c[2] ) + creal( c[2] ) * cimag( b[0] ) ) \
       + creal( b[2] ) * creal( c[3] ) - cimag( b[2] ) * cimag( c[3] ) - I * ( creal( b[2] ) * cimag( c[3] ) + creal( c[3] ) * cimag( b[2] ) ) ; 

  //a[2] = b[1]*conj( c[0] ) + b[3]*conj( c[1] )
  a[2] = -conj( a[1] ) ;  
  // a[3] = b[1]*conj( c[2] ) + b[3]*conj*( c[3] )
  a[3] = conj( a[0] ) ; 
#else
  return multab_dagdag( a , b , c ) ;
#endif
  return ;
}
#endif

#endif // immintrin.h
