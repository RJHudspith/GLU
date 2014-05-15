/*
    Copyright 2013 Renwick James Hudspith

    This file (solver.h) is part of GLU.

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
   @file solver.h
   @brief function protypes used for calculating eigenvalues of small matrices
 */
#ifndef GLU_SOLVER_H
#define GLU_SOLVER_H

/**
   @fn void squarert( double complex *__restrict res , const double complex z , const double complex R )
   @brief complex square root using de-Moivre's
   @param z :: the complex number being rooted
   @param R :: some parameter 
   @param res :: result of the rooting

   CALCULATION OF THE COMPLEX SQUARE ROOT WITH CONDITION 
   \f[ 
   \Re \left( R*\sqrt{ R^2 + g^3} \right) >= 0 
   \f]
   <br>
   very slow apparently 
 */
void 
squarert( double complex *__restrict res , 
	  const double complex z ,
	  const double complex R ) ;

/**
   @fn void cubert( double complex *__restrict res , const double complex z )
   @brief computes the cube root using de-Moivre's
   @param z :: complex number being rooted
   @param res :: result passed by reference
 */
void 
cubert( double complex *__restrict res ,
	const double complex z ) ;

/**
   @fn void Eigenvalues( GLU_complex z[ NC ] , const GLU_complex *__restrict U )
   @brief general eigenvalue calculator fo NCxNC small matrices
   @param U :: some complex square matrix
   @param z :: its eigenvalues
 */
void 
Eigenvalues( double complex z[ NC ] ,
	     const GLU_complex *__restrict U ) ;

/**
   @fn void Eigenvalues_suNC( double complex z[ NC ] , const GLU_complex U[ NCNC ] )
   @brief general eigenvalue calculator fo SU(NC) small matrices
   @param U :: some link matrix
   @param z :: its eigenvalues
   @warning U must be special unitary
   @return a flag saying whether we are close to eigenvlue degeneracy
 */
int 
Eigenvalues_suNC( double complex z[ NC ] ,
		 const GLU_complex U[ NCNC ] ) ;


/**
   @fn void Eigenvalues_hermitian( double z[ NC ] , const GLU_complex U[ NCNC ] )
   @brief general eigenvalue calculator for hermitian small matrices
   @param U :: some link matrix
   @param z :: its eigenvalues
   @warning U must be hermitian
 */
void 
Eigenvalues_hermitian( double z[ NC ] ,
		       const GLU_complex U[ NCNC ] ) ;

#endif
