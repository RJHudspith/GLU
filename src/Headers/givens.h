/*
    Copyright 2013 Renwick James Hudspith

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
   @fn void compute_pertinent_indices( void )
   @brief computes the various su2 subgoup indices of our NCxNC matrix
   @warning allocates memory

   I would like to thank Chroma, and in particular Edwards for this code.
 */
void
compute_pertinent_indices( void ) ;

/**
   @fn void free_su2_data( void )
   @brief frees the su2 subgroup struct allocated in givens.c
 */
void
free_su2_data( void ) ;

/**
   @fn void givens_reunit( double complex U[ NCNC ] ) ;
   @brief reunitarises a matrix with trace maximisation
   @param U :: matrix being reunitarised

   @warning overwrites the matrix U
 */
void
givens_reunit( GLU_complex U[ NCNC ] ) ;

#endif
