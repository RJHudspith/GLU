/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (MMULdag.c) is part of GLU.

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
   @file MMULdag.c
   @brief computes \f$ a = b^{\dagger}\times c \f$
 */
#include "Mainfile.h"

#if !( defined HAVE_IMMINTRIN_H ) || ( defined SINGLE_PREC )

#ifndef multabdag
// 3x3 mult a=b^{\dagger}.c//
void 
multabdag( GLU_complex a[ NCNC ] , 
	   const GLU_complex b[ NCNC ] , 
	   const GLU_complex c[ NCNC ] )
{
#if NC==3
  a[0] = conj( b[0] ) * c[0] + conj( b[3] ) * c[3] + conj( b[6] ) * c[6] ; \
  a[1] = conj( b[0] ) * c[1] + conj( b[3] ) * c[4] + conj( b[6] ) * c[7] ; \
  a[2] = conj( b[0] ) * c[2] + conj( b[3] ) * c[5] + conj( b[6] ) * c[8] ; \
  a[3] = conj( b[1] ) * c[0] + conj( b[4] ) * c[3] + conj( b[7] ) * c[6] ; \
  a[4] = conj( b[1] ) * c[1] + conj( b[4] ) * c[4] + conj( b[7] ) * c[7] ; \
  a[5] = conj( b[1] ) * c[2] + conj( b[4] ) * c[5] + conj( b[7] ) * c[8] ; \
  a[6] = conj( b[2] ) * c[0] + conj( b[5] ) * c[3] + conj( b[8] ) * c[6] ; \
  a[7] = conj( b[2] ) * c[1] + conj( b[5] ) * c[4] + conj( b[8] ) * c[7] ; \
  a[8] = conj( b[2] ) * c[2] + conj( b[5] ) * c[5] + conj( b[8] ) * c[8] ;
#elif NC==2
  a[0] = conj( b[0] ) * c[0] + conj( b[2] ) * c[2] ;	\
  a[1] = conj( b[0] ) * c[1] + conj( b[2] ) * c[3] ;	\
  a[2] = conj( b[1] ) * c[0] + conj( b[3] ) * c[2] ;	\
  a[3] = conj( b[1] ) * c[1] + conj( b[3] ) * c[3] ;
#else
  size_t i , j , m ;
  register GLU_complex sum ;
  register GLU_real REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = 0.0 ;
      for( m = 0 ; m < NC ; m++ ) {
	REB = creal( b[ i + NC*m ] ) ; IMB = cimag( b[ i + NC*m ] ) ;
	REC = creal( c[ j + NC*m ] ) ; IMC = cimag( c[ j + NC*m ] ) ;
	sum += REB * REC + IMB * IMC + I * ( REB * IMC - IMB * REC ) ;
      }
      a[ j + NC*i ] = sum ;
    }
  }
#endif
  return ;
}
#endif // !defined multabdag

#endif // immintrin.h
