/*
    Copyright 2013 Renwick James Hudspith

    This file (MMULdagdag.h) is part of GLU.

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
   @file MMULdagdag.h
   @brief prototype functions for matrix multiplication of generic matrices and SU(NC) only
 */
#ifndef GLU_MMULDAGDAG_H
#define GLU_MMULDAGDAG_H

#ifndef GLU_BGQ

/**
   @fn void multab_dagdag( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief computes \f$ a = b^{\dagger} \times c^{\dagger} \f$

   loop unrolled, standard matrix multiply
 **/
void 
multab_dagdag( GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC ] , 
	       const GLU_complex c[ NCNC ] ) ;

#if NC < 4

/**
   @fn void multab_dagdag_suNC( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief computes \f$ a = b^{\dagger} \times c^{\dagger}\f$ sped up for suNC

   This relies on the fact that computing the bottom row's minor's requires
   less operations than computing the full matrix multiply also uses hand unrolled complex-complex multiply

   @warning only use-able for \f$ a,b,c \in SU(NC) \f$
 **/
void 
multab_dagdag_suNC( GLU_complex a[ NCNC ] , 
		    const GLU_complex b[ NCNC ] , 
		    const GLU_complex c[ NCNC ] ) ;

#else
  #define multab_dagdag_suNC multab_dagdag
#endif

#else // GLU_BGQ ifdef

#if NC==3
#define multab_dagdag( a , b , c )					\
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
#define multab_dagdag( a , b , c )					\
  a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[2] ) * conj( c[1] )  ;	\
  a[1] = conj( b[0] ) * conj( c[2] ) + conj( b[2] ) * conj( c[3] )  ;	\
  a[2] = conj( b[1] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] )  ;	\
  a[3] = conj( b[1] ) * conj( c[2] ) + conj( b[3] ) * conj( c[3] )  ; 
#else
 void 
 multab_dagdag( GLU_complex a[ NCNC ] , 
		const GLU_complex b[ NCNC ] , 
		const GLU_complex c[ NCNC ] ) ;
#endif

// SUNC variants
#if NC==3
#define multab_dagdag_suNC( a , b , c )		\
  a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] ) + conj( b[6] ) * conj( c[2] ) ; \
  a[1] = conj( b[0] ) * conj( c[3] ) + conj( b[3] ) * conj( c[4] ) + conj( b[6] ) * conj( c[5] ) ; \
  a[2] = conj( b[0] ) * conj( c[6] ) + conj( b[3] ) * conj( c[7] ) + conj( b[6] ) * conj( c[8] ) ; \
  a[3] = conj( b[1] ) * conj( c[0] ) + conj( b[4] ) * conj( c[1] ) + conj( b[7] ) * conj( c[2] ) ; \
  a[4] = conj( b[1] ) * conj( c[3] ) + conj( b[4] ) * conj( c[4] ) + conj( b[7] ) * conj( c[5] ) ; \
  a[5] = conj( b[1] ) * conj( c[6] ) + conj( b[4] ) * conj( c[7] ) + conj( b[7] ) * conj( c[8] ) ; \
  a[6] = conj( a[1] * a[5] - a[2] * a[4] ) ;			\
  a[7] = conj( a[2] * a[3] - a[0] * a[5] ) ;			\
  a[8] = conj( a[0] * a[4] - a[1] * a[3] ) ; 
#elif NC==2
#define multab_dagdag_suNC( a , b , c )		\
  a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[2] ) * conj( c[1] ) ; \
  a[1] = conj( b[0] ) * conj( c[2] ) + conj( b[2] ) * conj( c[3] ) ; \
  a[2] = -conj( a[1] ) ;					     \
  a[3] = conj( a[0] ) ; 
#else
  #define multab_dagdag_suNC multab_dagdag
#endif

#endif // end of the BGQ inlining loop ... 

#endif
