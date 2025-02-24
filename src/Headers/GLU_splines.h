/*
Copyright 2013-2025 Renwick James Hudspith

    This file (GLU_splines.h) is part of GLU.

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
   @file GLU_splines.h
   @brief prototype functions for cubic spline evaluations
 */
#ifndef GLU_SPLINES_H
#define GLU_SPLINES_H

/**
   @fn double cubic_eval( const double *__restrict x , const double *__restrict y , const double *__restrict der , const double mu , const size_t datalength )
   @brief computes the cubic spline evaluation at the point "mu"
   @param x :: x-data
   @param y :: y-data
   @param der :: the derivatives dy/dx for each x point
   @param mu :: interpolation point
   @param datalength :: the length of the y,x and der arrays
 */
double
cubic_eval( const double *__restrict x ,
	    const double *__restrict y ,
	    const double *__restrict der ,
	    const double mu ,
	    const size_t datalength ) ;

/**
   @fn double cubic_min( const double *__restrict x , const double *__restrict y , const double *__restrict der , const size_t change_up )
   @brief evaluates the minimum of a cubic spline interpolation
   @param x :: x-data
   @param y :: y-data
   @param der :: the derivatives dy/dx for each x point
   @param change up :: the index of y for which the derivative changes sign

   @warning requires the minimum to be bounded and have change_up as the upper index for where the derivative changes
 */
double
cubic_min( const double *__restrict x ,
	   const double *__restrict y ,
	   const double *__restrict der ,
	   const size_t change_up ) ;

/**
   @fn void spline_derivative( double *__restrict der , const double *__restrict x , const double *__restrict y , const int N )
   @brief computes dy/dx at each x-point
   @param der :: the derivative; must have workspace of length N
   @param x :: the x-data
   @param y :: y-data
   @param N :: length of the x and y data

   @warning If N is less than 5 we use lower (than 4th) order finite difference defs
 */
void
spline_derivative( double *__restrict der ,
		   const double *__restrict x ,
		   const double *__restrict y ,
		   const size_t N ) ;

/**
   @fn double solve_spline( const double *__restrict x , const double *__restrict y , const double *__restrict der , const double mu , const size_t datalength )
   @brief solve the cubic spline y(x) = mu
 */
double
solve_derspline( const double *__restrict x ,
		 const double *__restrict y ,
		 const double *__restrict der ,
		 const double mu ,
		 const size_t change_up ) ;

/**
   @fn double solve_spline( const double *__restrict x , const double *__restrict y , const double *__restrict der , const double mu , const size_t datalength )
   @brief solve the cubic spline y(x) = mu
 */
double
solve_spline( const double *__restrict x ,
	      const double *__restrict y ,
	      const double *__restrict der ,
	      const double mu ,
	      const size_t change_up ) ;

#endif
