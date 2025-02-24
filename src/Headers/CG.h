/*
Copyright 2013-2025 Renwick James Hudspith

    This file (CG.h) is part of GLU.

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
   @file CG.h
   @brief common routines between Landau and Coulomb gauge fixing
 */
#ifndef GLU_CG_H
#define GLU_CG_H

/**
   @fn double PRfmax( const double a , const double b )
   @brief returns the maximum of a and b
 */
double
PRfmax( const double a , const double b ) ;

/**
   @fn void set_gauge_matrix( GLU_complex *__restrict gx , const GLU_complex *__restrict *__restrict in , const double alpha , const size_t i ) 
   @brief unpacks and exponentiates the derivative into the array gx at site i
   @param gx :: gauge transformation matrix
   @param in :: derivative at site i
   @param i :: site index
 */
void
set_gauge_matrix( GLU_complex *__restrict gx ,
		  const GLU_complex *__restrict *__restrict in ,
		  const double alpha ,
		  const size_t i ) ;

#endif
