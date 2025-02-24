/*
Copyright 2013-2025 Renwick James Hudspith

    This file (invert.h) is part of GLU.

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
   @file invert.h
   @brief function definition for taking the naive numerical inverse of a matrix
   Uses Gauss-Jordan with pivoting unless specified otherwise or lapacke.h is used
*/
#ifndef GLU_INVERT_H
#define GLU_INVERT_H

/**
   @fn int inverse( GLU_complex M_1[ NCNC ] , const GLU_complex M[ NCNC ] )
   @brief inverts a matrix using the usual minors definition dividing the determinant
   @param M_1 :: the inverse of the matrix "M"
   @param M :: the matrix being inverted 

   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int 
inverse( GLU_complex M_1[ NCNC ] , 
	 const GLU_complex M[ NCNC ] ) ;

/**
   @fn void newton_approx_inverse( GLU_complex Zinv[ NCNC ] , const GLU_complex Z[ NCNC ] )
   @brief computes a newton iteration for the inverse, it is not very good
   @param Zinv :: the matrix inverse
   @param Z :: the matrix that is having its inverse taken
 */
void
newton_approx_inverse( GLU_complex Zinv[ NCNC ] ,
		       const GLU_complex Z[ NCNC ] ) ;

#endif
