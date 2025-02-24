/*
Copyright 2013-2025 Renwick James Hudspith

    This file (givens.h) is part of GLU.

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
   @file givens.h
   @brief prototype functions for the trace maximisation using givens rotations

   some credit should go to Urs Wenger and ETMC and Edwards and Chroma
 */
#ifndef GLU_GIVENS_H
#define GLU_GIVENS_H

/**
   @fn void givens_reunit( double complex U[ NCNC ] ) ;
   @brief reunitarises a matrix with trace maximisation
   @param U :: matrix being reunitarised

   @warning overwrites the matrix U
 */
void
givens_reunit( GLU_complex U[ NCNC ] ) ;

/**
   @fn void OrRotation( GLU_complex *s0 , GLU_complex *s1 , const GLU_complex U[ NCNC ] , const double OrParam , const int su2_index )
   @brief compute the two defining parameters of the su(2) representation
   @param U :: matrix whose subgroup is computed
   @param s0 :: top left su(2) element
   @param s1 :: top right su(2) element
   @param OrParam :: overrelaxation parameter
   @param su2_index :: su(2) subgroup index
 */
void
OrRotation( GLU_complex *s0 , 
	    GLU_complex *s1 ,
	    const GLU_complex U[ NCNC ] , 
	    const double OrParam ,
	    const size_t su2_index ) ;

#endif
