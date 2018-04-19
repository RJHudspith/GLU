/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (MMUL_SUNC.h) is part of GLU.

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
   @file MMUL_SUNC.h
   @brief prototype functions for matrix multiplication routines
 */
#ifndef GLU_MMUL_SUNC_H
#define GLU_MMUL_SUNC_H

#ifndef GLU_BGQ // turns on the inline matrix multiplies

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
  #include "BGQ_mmuls.h"
#endif // end of the inline matrix mul loop
#endif 
