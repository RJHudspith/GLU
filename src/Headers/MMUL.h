/*
Copyright 2013-2025 Renwick James Hudspith

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

#else // else we inline the matrix multiplies ...
  #include "BGQ_mmuls.h"
#endif

#endif 
