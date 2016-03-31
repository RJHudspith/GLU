/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (vandermonde.h) is part of GLU.

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
   @file vandermonde.h
   @brief prototype function for solving a complex vandermonde system with real eigenvalues
 */
#ifndef GLU_VANDERMONDE_H
#define GLU_VANDERMONDE_H

/**
   \fn void solvandermonde( double complex f[ NC ] , const double x[ NC ] )
   \brief solves the vandermonde system using interpolating polynomials
   @param f :: the f constants
   @param x :: the real eigenvalues
   @warning numerically pretty unstable, f's overwritten
   see e.g Gloub and van Loan <a href="http://books.google.co.uk/books/about/Matrix_Computations.html?id=mlOa7wPX6OYC&redir_esc=y" > book </a>
 */
void 
solvandermonde( double complex f[ NC ] , 
		const double x[ NC ] ) ;


/**
   \fn void solvandermonde_new( const double q[ NC ] , double complex f[ NC ] )
   \brief solves the vandermonde system using interpolating polynomials
   @param q :: the real eigenvalues
   @param f :: the f constants
   @warning numerically quasi-unstable, f's overwritten
*/
void 
solvandermonde_new( const double q[ NC ] , 
		    double complex f[ NC ] ) ;
#endif
