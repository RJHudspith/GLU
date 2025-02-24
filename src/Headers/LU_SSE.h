/*
Copyright 2013-2025 Renwick James Hudspith

    This file (LU_SSE.h) is part of GLU.

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
   @file LU_SSE.h
   @brief LU decomposition determinant of an complex matrix (SSE2 variant)
 */
#ifndef GLU_LU_SSE_H
#define GLU_LU_SSE_H

#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC ) && (NC>3)

/**
   @fn double complex LU_det( const size_t N , const GLU_complex U[ N*N ] )
   @brief computes the determinant of a matrix U
   @param N :: Side length of the square matrix
   @param U :: Matrix having its determinant taken
   @return the determinant
 */
double complex
LU_det( const size_t N , 
	const GLU_complex U[ N*N ] ) ;

/**
   @fn double complex LU_det_overwrite( const size_t N , GLU_complex U[ N*N ] )
   @brief computes the determinant of a matrix U
   @param N :: Side length of the square matrix
   @param U :: Matrix having its determinant taken, overwritten
   @return the determinant
   @warning overwrites space in U
 */
double complex
LU_det_overwrite( const size_t N , 
		  GLU_complex U[ N*N ] ) ;

#endif

#endif
