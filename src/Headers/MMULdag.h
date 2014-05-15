/*
    Copyright 2013 Renwick James Hudspith

    This file (MMULdag.h) is part of GLU.

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
   @file MMULdag.h
   @brief prototype functions for multiplication for \f$ a = b^{\dagger}\times c\f$ 
 */

#ifndef GLU_MMULDAG_H
#define GLU_MMULDAG_H

#ifndef GLU_BGQ // dirty, force inlining for the Q

/**
   @fn void multabdag( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief Computes general Nc x Nc matrix multiplication with \f$ a = b^{\dagger}\times c \f$

   general and slow, uses the complex multiplication defined in complex.h
 **/
void 
multabdag( GLU_complex a[ NCNC ] , 
	   const GLU_complex b[ NCNC ] , 
	   const GLU_complex c[ NCNC ] ) ;

#if NC < 4

/**
   @fn void multabdag_suNC( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief Computes general Nc x Nc matrix multiplication with \f$ a = b^{\dagger}\times c \f$

   sped up for suNC by computing the signed minors for the bottom row and hand-expanded complex-complex multiplications

   @warning only use-able for \f$ a,b,c \in SU(NC) \f$
 **/
void 
multabdag_suNC( GLU_complex a[ NCNC ] , 
		const GLU_complex b[ NCNC ] , 
		const GLU_complex c[ NCNC ] ) ;

#else
  #define multabdag_suNC multabdag
#endif

#else

#if NC==3
    #define multabdag( a , b , c )					     \
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
  #define multabdag( a , b , c )				\
    a[0] = conj( b[0] ) * c[0] + conj( b[2] ) * c[2] ;	\
    a[1] = conj( b[0] ) * c[1] + conj( b[2] ) * c[3] ;	\
    a[2] = conj( b[1] ) * c[0] + conj( b[3] ) * c[2] ;	\
    a[3] = conj( b[1] ) * c[1] + conj( b[3] ) * c[3] ;
  #else
  void 
  multabdag( GLU_complex a[ NCNC ] ,
	     const GLU_complex b[ NCNC ] , 
	     const GLU_complex c[ NCNC ] ) ;
  #endif

  #if NC==3
    #define multabdag_suNC( a , b , c )					\
    a[0] = conj( b[0] ) * c[0] + conj( b[3] ) * c[3] + conj( b[6] ) * c[6] ; \
    a[1] = conj( b[0] ) * c[1] + conj( b[3] ) * c[4] + conj( b[6] ) * c[7] ; \
    a[2] = conj( b[0] ) * c[2] + conj( b[3] ) * c[5] + conj( b[6] ) * c[8] ; \
    a[3] = conj( b[1] ) * c[0] + conj( b[4] ) * c[3] + conj( b[7] ) * c[6] ; \
    a[4] = conj( b[1] ) * c[1] + conj( b[4] ) * c[4] + conj( b[7] ) * c[7] ; \
    a[5] = conj( b[1] ) * c[2] + conj( b[4] ) * c[5] + conj( b[7] ) * c[8] ; \
    a[6] = conj( a[1]  *  a[5] - a[2]  *  a[4] ) ;			\
    a[7] = conj( a[2]  *  a[3] - a[0]  *  a[5] ) ;			\
    a[8] = conj( a[0]  *  a[4] - a[1]  *  a[3] ) ;
  #elif NC==2
    #define multabdag_suNC( a , b , c )		\
    a[0] = conj( b[0] ) * c[0] + conj( b[2] ) * c[2] ;	\
    a[1] = conj( b[0] ) * c[1] + conj( b[2] ) * c[3] ;	\
    a[2] = -conj( a[1] ) ;				\
    a[3] = conj( a[0] ) ;
  #else
    #define multabdag_suNC multabdag
  #endif

#endif

#endif
