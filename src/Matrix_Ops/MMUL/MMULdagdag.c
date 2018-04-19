/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (MMULdagdag.c) is part of GLU.

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
   @file MMULdagdag.c
   @brief computes \f$ a = b^{\dagger} c^{\dagger} \f$ matrix multiplication
 */
#include "Mainfile.h"

#if !( defined HAVE_IMMINTRIN_H ) || ( defined SINGLE_PREC )

#ifndef multab_dagdag
// a = b^{\dagger} * c^{\dagger}
void 
multab_dagdag( GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC ] , 
	       const GLU_complex c[ NCNC ] )
{
#if NC==3
  a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] ) + conj( b[6] ) * conj( c[2] ) ; \
  a[1] = conj( b[0] ) * conj( c[3] ) + conj( b[3] ) * conj( c[4] ) + conj( b[6] ) * conj( c[5] ) ; \
  a[2] = conj( b[0] ) * conj( c[6] ) + conj( b[3] ) * conj( c[7] ) + conj( b[6] ) * conj( c[8] ) ; \
  a[3] = conj( b[1] ) * conj( c[0] ) + conj( b[4] ) * conj( c[1] ) + conj( b[7] ) * conj( c[2] ) ; \
  a[4] = conj( b[1] ) * conj( c[3] ) + conj( b[4] ) * conj( c[4] ) + conj( b[7] ) * conj( c[5] ) ; \
  a[5] = conj( b[1] ) * conj( c[6] ) + conj( b[4] ) * conj( c[7] ) + conj( b[7] ) * conj( c[8] ) ; \
  a[6] = conj( b[2] ) * conj( c[0] ) + conj( b[5] ) * conj( c[1] ) + conj( b[8] ) * conj( c[2] ) ; \
  a[7] = conj( b[2] ) * conj( c[3] ) + conj( b[5] ) * conj( c[4] ) + conj( b[8] ) * conj( c[5] ) ; \
  a[8] = conj( b[2] ) * conj( c[6] ) + conj( b[5] ) * conj( c[7] ) + conj( b[8] ) * conj( c[8] ) ; 
#elif NC==2
  a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[2] ) * conj( c[1] )  ;	\
  a[1] = conj( b[0] ) * conj( c[2] ) + conj( b[2] ) * conj( c[3] )  ;	\
  a[2] = conj( b[1] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] )  ;	\
  a[3] = conj( b[1] ) * conj( c[2] ) + conj( b[3] ) * conj( c[3] )  ; 
#else
  size_t i , j , m ;
  register GLU_complex sumr = 0.0 , sumi = 0.0 ;
  register GLU_real REB , IMB , REC , IMC ;
  const GLU_complex *pC , *pB ;
  for( i = 0 ; i < NC ; i++ ) {
    pC = c ;
    for( j = 0 ; j < NC ; j++ ) {
      pB = b ;
      sumr = sumi = 0.0 ;
      for( m = 0 ; m < NC ; m++ ) {
	REB = creal( pB[i] ) ; IMB = cimag( pB[i] ) ;
	REC = creal( pC[m] ) ; IMC = cimag( pC[m] ) ;
	sumr += REB * REC - IMB * IMC ;
	sumi += REB * IMC + IMB * REC ;
	pB += NC ;
      }
      a[ j + NC*i ] = sumr + I * sumi ;
      pC += NC ;
    }
  }
#endif
  return ;
}
#endif // !defined multabdagdag

#endif // immintrin.h
