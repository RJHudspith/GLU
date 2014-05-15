/*
    Copyright 2013 Renwick James Hudspith

    This file (MMUL.h) is part of GLU.

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
   @file MMUL.h
   @brief prototype functions for matrix multiplication routines
 */

#ifndef GLU_MMUL_H
#define GLU_MMUL_H

/**
   @fn void multab_atomic_left( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] )
   @brief NC x NC matrix multiply a = b * a

   computes \f$ a = b \times a \f$ general NC x NC matrix multiply does not require specific types of matrices
   @warning overwrites the matix @a
 **/
void 
multab_atomic_left( GLU_complex a[ NCNC ] , 
		    const GLU_complex b[ NCNC ] ) ;

/**
   @fn void multab_atomic_right( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] )
   @brief NC x NC matrix multiply a = a * b

   computes \f$ a = a \times b \f$ general NC x NC matrix multiply does not require specific types of matrices
   @warning overwrites the matix @a
 **/
void 
multab_atomic_right( GLU_complex a[ NCNC ] , 
		     const GLU_complex b[ NCNC ] ) ;

#ifndef GLU_BGQ // turns on the inline matrix multiplies

/**
   @fn void multab( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief NC x NC matrix multiply

   computes \f$ a = b \times c \f$ general Nc x Nc matrix multiply does not require specific types of matrices
 **/
void 
multab( GLU_complex a[ NCNC ] , 
	const GLU_complex b[ NCNC ] , 
	const GLU_complex c[ NCNC ] ) ;

#if NC < 4
/**
   @fn void multab_suNC( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief SU(Nc)- tuned matrix multiply 

   computes \f$ a = b \times c \f$  tuned for suNC by computing the signed minors for the bottom row

   @warning only use-able for \f$ a,b,c \in SU(NC) \f$
 **/
void 
multab_suNC( GLU_complex a[ NCNC ] , 
	     const GLU_complex b[ NCNC ] , 
	     const GLU_complex c[ NCNC ] ) ;
#else
  #define multab_suNC multab
#endif

#else // else we inline the matrix multiplies ...

  // This was not helpful! -> Perhaps bgq?
  #if NC==3
    #define multab( a , b , c )				\
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
    #define multab( a , b , c )			\
    a[0] = b[0] * c[0] + b[1] * c[2] ;		\
    a[1] = b[0] * c[1] + b[1] * c[3] ;		\
    a[2] = b[2] * c[0] + b[3] * c[2] ;		\
    a[3] = b[2] * c[1] + b[3] * c[3] ;		 
  #else
    // matrix function call because the naive routine is expensive anyway 
    void 
    multab( GLU_complex a[ NCNC ] , 
	    const GLU_complex b[ NCNC ] , 
	    const GLU_complex c[ NCNC ] ) ;
  #endif

  #if NC==3 
    #define multab_suNC( a , b , c )			\
    a[0] = b[0] * c[0] + b[1] * c[3] + b[2] * c[6] ;	\
    a[1] = b[0] * c[1] + b[1] * c[4] + b[2] * c[7] ;	\
    a[2] = b[0] * c[2] + b[1] * c[5] + b[2] * c[8] ;	\
    a[3] = b[3] * c[0] + b[4] * c[3] + b[5] * c[6] ;	\
    a[4] = b[3] * c[1] + b[4] * c[4] + b[5] * c[7] ;	\
    a[5] = b[3] * c[2] + b[4] * c[5] + b[5] * c[8] ;	\
    a[6] = conj( a[1] * a[5] - a[2] * a[4] ) ;	\
    a[7] = conj( a[2] * a[3] - a[0] * a[5] ) ;	\
    a[8] = conj( a[0] * a[4] - a[1] * a[3] ) ;	
  #elif NC==2
    #define multab_suNC( a , b , c )		\
    a[0] = b[0] * c[0] + b[1] * c[2] ;		\
    a[1] = b[0] * c[1] + b[1] * c[3] ;		\
    a[2] = -conj( a[1] ) ;			\
    a[3] = conj( a[0] ) ;
  #else
    #define multab_suNC multab
  #endif

#endif // end of the inline matrix mul loop

#endif 
