/*
    Copyright 2013 Renwick James Hudspith

    This file (MMUL_dag.c) is part of GLU.

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
   @file MMUL_dag.c
   @brief matrix multiply \f$ a = b \times c^{\dagger}\f$ , where \f$ a,b,c \in NCxNC \f$ GLU_complex matrices 
 */

#include "Mainfile.h"

#if !( defined HAVE_IMMINTRIN_H ) || ( defined SINGLE_PREC )

// a = b * c^{\dagger}
#ifndef multab_dag
void 
multab_dag( GLU_complex a[ NCNC ] , 
	    const GLU_complex b[ NCNC ] , 
	    const GLU_complex c[ NCNC ] )
{
#if NC==3
  a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) + b[2] * conj( c[2] ) ; \
  a[1] = b[0] * conj( c[3] ) + b[1] * conj( c[4] ) + b[2] * conj( c[5] ) ; \
  a[2] = b[0] * conj( c[6] ) + b[1] * conj( c[7] ) + b[2] * conj( c[8] ) ; \
  a[3] = b[3] * conj( c[0] ) + b[4] * conj( c[1] ) + b[5] * conj( c[2] ) ; \
  a[4] = b[3] * conj( c[3] ) + b[4] * conj( c[4] ) + b[5] * conj( c[5] ) ; \
  a[5] = b[3] * conj( c[6] ) + b[4] * conj( c[7] ) + b[5] * conj( c[8] ) ; \
  a[6] = b[6] * conj( c[0] ) + b[7] * conj( c[1] ) + b[8] * conj( c[2] ) ; \
  a[7] = b[6] * conj( c[3] ) + b[7] * conj( c[4] ) + b[8] * conj( c[5] ) ; \
  a[8] = b[6] * conj( c[6] ) + b[7] * conj( c[7] ) + b[8] * conj( c[8] ) ; 
#elif NC==2
  a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) ;	\
  a[1] = b[0] * conj( c[2] ) + b[1] * conj( c[3] ) ;	\
  a[2] = b[2] * conj( c[0] ) + b[3] * conj( c[1] ) ;	\
  a[3] = b[2] * conj( c[2] ) + b[3] * conj( c[3] ) ;
#else // instead of inlining we have a function call
  #if NC%2==0
  int i , j , m ;
  const GLU_complex *pB = b , *pC ;
  register GLU_complex sum , sum2 ;
  register GLU_real REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    pC = c ;
    for( j = 0 ; j < NC ; j+=2 ) {
      sum = sum2 = 0. ;
      for( m = 0 ; m < NC ; m++ ) {
	// set REB
	REB = creal( pB[m] ) ; IMB = cimag( pB[m] ) ;
	REC = creal( pC[m] ) ; IMC = cimag( pC[m] ) ;
	sum += REB * REC + IMB * IMC + I * ( REC * IMB - REB * IMC ) ;
	REC = creal( pC[m+NC] ) ; IMC = cimag( pC[m+NC] ) ;
	sum2 += REB * REC + IMB * IMC + I * ( REC * IMB - REB * IMC ) ;
      }
      a[ j + NC*i ] = sum ;
      a[ j + 1 + NC*i ] = sum2 ;
      // double increment
      pC += NC ;
      pC += NC ;
    }
    pB += NC ; 
  }
  #else
  int i , j , m ;
  register GLU_complex sum ;
  register GLU_real REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = 0. ;
      for( m = 0 ; m < NC ; m++ ) {
	REB = creal( b[ m + NC*i ] ) ; IMB = cimag( b[ m + NC*i ] ) ;
	REC = creal( c[ m + NC*j ] ) ; IMC = cimag( c[ m + NC*j ] ) ;
	sum += REB * REC + IMB * IMC + I * ( REC * IMB - REB * IMC ) ;
      }
      a[ j + NC*i ] = sum ;
    }
  }
  #endif
#endif
  return ;
}
#endif

#ifndef multab_dag_suNC

void //__attribute__((hot))
multab_dag_suNC( GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] , 
		 const GLU_complex c[ NCNC ] )
{
#if NC == 3
  //a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) + b[2] * conj( c[2] ) ; 
  a[0] = creal( b[0] ) * creal( c[0] ) + cimag( b[0] ) * cimag( c[0] ) - I * ( creal( b[0] ) * cimag( c[0] ) - creal( c[0] ) * cimag( b[0] ) ) \
       + creal( b[1] ) * creal( c[1] ) + cimag( b[1] ) * cimag( c[1] ) - I * ( creal( b[1] ) * cimag( c[1] ) - creal( c[1] ) * cimag( b[1] ) ) \
       + creal( b[2] ) * creal( c[2] ) + cimag( b[2] ) * cimag( c[2] ) - I * ( creal( b[2] ) * cimag( c[2] ) - creal( c[2] ) * cimag( b[2] ) ) ; 

  //a[1] = b[0] * conj( c[3] ) + b[1] * conj( c[4] ) + b[2] * conj( c[5] ) ; 
  a[1] = creal( b[0] ) * creal( c[3] ) + cimag( b[0] ) * cimag( c[3] ) - I * ( creal( b[0] ) * cimag( c[3] ) - creal( c[3] ) * cimag( b[0] ) ) \
       + creal( b[1] ) * creal( c[4] ) + cimag( b[1] ) * cimag( c[4] ) - I * ( creal( b[1] ) * cimag( c[4] ) - creal( c[4] ) * cimag( b[1] ) ) \
       + creal( b[2] ) * creal( c[5] ) + cimag( b[2] ) * cimag( c[5] ) - I * ( creal( b[2] ) * cimag( c[5] ) - creal( c[5] ) * cimag( b[2] ) ) ; 

  //a[2] = b[0] * conj( c[6] ) + b[1] * conj( c[7] ) + b[2] * conj( c[8] ) ; 
  a[2] = creal( b[0] ) * creal( c[6] ) + cimag( b[0] ) * cimag( c[6] ) - I * ( creal( b[0] ) * cimag( c[6] ) - creal( c[6] ) * cimag( b[0] ) ) \
       + creal( b[1] ) * creal( c[7] ) + cimag( b[1] ) * cimag( c[7] ) - I * ( creal( b[1] ) * cimag( c[7] ) - creal( c[7] ) * cimag( b[1] ) ) \
       + creal( b[2] ) * creal( c[8] ) + cimag( b[2] ) * cimag( c[8] ) - I * ( creal( b[2] ) * cimag( c[8] ) - creal( c[8] ) * cimag( b[2] ) ) ; 

  // middle row //
  //a[3] = b[3] * conj( c[0] ) + b[4] * conj( c[1] ) + b[5] * conj( c[2] ) ; 
  a[3] = creal( b[3] ) * creal( c[0] ) + cimag( b[3] ) * cimag( c[0] ) - I * ( creal( b[3] ) * cimag( c[0] ) - creal( c[0] ) * cimag( b[3] ) ) \
       + creal( b[4] ) * creal( c[1] ) + cimag( b[4] ) * cimag( c[1] ) - I * ( creal( b[4] ) * cimag( c[1] ) - creal( c[1] ) * cimag( b[4] ) ) \
       + creal( b[5] ) * creal( c[2] ) + cimag( b[5] ) * cimag( c[2] ) - I * ( creal( b[5] ) * cimag( c[2] ) - creal( c[2] ) * cimag( b[5] ) ) ; 

  //a[4] = b[3] * conj( c[3] ) + b[4] * conj( c[4] ) + b[5] * conj( c[5] ) ; 
  a[4] = creal( b[3] ) * creal( c[3] ) + cimag( b[3] ) * cimag( c[3] ) - I * ( creal( b[3] ) * cimag( c[3] ) - creal( c[3] ) * cimag( b[3] ) ) \
       + creal( b[4] ) * creal( c[4] ) + cimag( b[4] ) * cimag( c[4] ) - I * ( creal( b[4] ) * cimag( c[4] ) - creal( c[4] ) * cimag( b[4] ) ) \
       + creal( b[5] ) * creal( c[5] ) + cimag( b[5] ) * cimag( c[5] ) - I * ( creal( b[5] ) * cimag( c[5] ) - creal( c[5] ) * cimag( b[5] ) ) ; 

  //a[5] = b[3] * conj( c[6] ) + b[4] * conj( c[7] ) + b[5] * conj( c[8] ) ;
  a[5] = creal( b[3] ) * creal( c[6] ) + cimag( b[3] ) * cimag( c[6] ) - I * ( creal( b[3] ) * cimag( c[6] ) - creal( c[6] ) * cimag( b[3] ) ) \
       + creal( b[4] ) * creal( c[7] ) + cimag( b[4] ) * cimag( c[7] ) - I * ( creal( b[4] ) * cimag( c[7] ) - creal( c[7] ) * cimag( b[4] ) ) \
       + creal( b[5] ) * creal( c[8] ) + cimag( b[5] ) * cimag( c[8] ) - I * ( creal( b[5] ) * cimag( c[8] ) - creal( c[8] ) * cimag( b[5] ) ) ; 

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
  //faster multab?   
  //a[0] = b[0]*conj( c[0] ) + b[1]*conj( c[2] )
  a[0] = creal( b[0] ) * creal( c[0] ) + cimag( b[0] ) * cimag( c[0] ) - I * ( creal( b[0] ) * cimag( c[0] ) - creal( c[0] ) * cimag( b[0] ) )   \
       + creal( b[1] ) * creal( c[1] ) + cimag( b[1] ) * cimag( c[1] ) - I * ( creal( b[1] ) * cimag( c[1] ) - creal( c[1] ) * cimag( b[1] ) ) ; \

  //a[1] = b[0] * conj( c[2] ) + b[1] * conj( c[3] ) 
  a[1] = creal( b[0] ) * creal( c[2] ) + cimag( b[0] ) * cimag( c[2] ) - I * ( creal( b[0] ) * cimag( c[2] ) - creal( c[2] ) * cimag( b[0] ) )   \
       + creal( b[1] ) * creal( c[3] ) + cimag( b[1] ) * cimag( c[3] ) - I * ( creal( b[1] ) * cimag( c[3] ) - creal( c[3] ) * cimag( b[1] ) ) ; \

  //a[2] = b[2]*conj( c[0] ) + b[3]*conj( c[1] )
  a[2] = -conj( a[1] ) ;  
  // a[3] = b[2]*conj( c[2] ) + b[3]*conj*( c[3] )
  a[3] = conj( a[0] ) ; 
#else
  return multab_dag( a , b , c ) ;
#endif
  return ;
}
#endif

#endif // immintrin.h
