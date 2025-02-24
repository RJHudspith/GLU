/*
Copyright 2013-2025 Renwick James Hudspith

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
   @file MMUL_triple.h
   @brief prototype functions for matrix multiplication routines
 */
#ifndef GLU_MMUL_TRIPLE_H
#define GLU_MMUL_TRIPLE_H

/**
   @fn void multabcdag_suNC( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] , const GLU_complex d[ NCNC ] )
   @brief SU(Nc) multiply a = b.c.d^\dagger
 */
void 
multabcdag_suNC( GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] , 
		 const GLU_complex c[ NCNC ] ,
		 const GLU_complex d[ NCNC ] ) ;


/**
   @fn void multadagbc_suNC( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] , const GLU_complex d[ NCNC ] )
   @brief SU(Nc) multiply a = b^\dagger.c.d
 */
void 
multadagbc_suNC( GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] , 
		 const GLU_complex c[ NCNC ] ,
		 const GLU_complex d[ NCNC ] ) ;

#endif
