/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (vandermonde.c) is part of GLU.

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
   @file vandermonde.c 
   @brief Iterative solution to the vandermonde system Va=f, which overwrites f with a
 */

#include "Mainfile.h"

// this code computes the solutions of a Vandermonde system overwritten in c[]//
// where the elements of the matrix are real and the solution f in Va=f are complex -> e^(iq_n)// 
void 
solvandermonde( double complex f[ NC ] ,
		const double x[ NC ] ) 
{
  /// compute the "newton representation" of interpolating polynomial ///
  size_t i , k ;
  for( k = 0 ; k < NC - 1 ; k++ ) { 
    for( i = NC - 1 ; i > k ; i-- ) {
      f[i] = ( f[i] - f[i-1] ) / ( x[i] - x[i-k-1] ) ;
    }
  }
  /// plug in values for the interpolator ///
  for( k = NC - 1 ; k != 0 ; k-- ) {
    for(i = k ; i < NC - 1 ; i++ ) {
      f[i] = f[i] - ( f[ i + 1 ] * x[k] ) ;
    }
  }
  return;
}

/*
  SOLVES THE SYSTEM

  | 1   q0   q0^2  ... q0^n | | f0 |   | e^{iq0} |
  | 1   q1   q1^2  ... q1^n | | f1 |   | e^{iq1} |
  | .   .     .    ...  .   | | .  | = | .       |
  | .   .     .    ...  .   | | .  | = | .       |
  | .   .     .    ...  .   | | .  |   | .       |
  | 1   qn   qn^2  ... qn^n | | fn |   | e^{iqn} |   

  Where fn is expected to be e^{iqn} when this function is called
 */
void 
solvandermonde_new( const double q[ NC ] , 
		    double complex f[ NC ] ) 
{
  /// compute the "newton representation" of interpolating polynomial ///
  size_t i , k ;
  for( k = 0 ; k < NC - 1 ; k++ ) { 
    for( i = NC - 1 ; i > k ; i-- ) {
      f[i] = ( f[i] - f[i-1] ) / ( q[i] - q[i-k-1] ) ;
    }
  }
  /// plug in values for the interpolator ///
  for( k = NC - 1 ; k != 0 ; k-- ) {
    for(i = k ; i < NC - 1 ; i++ ) {
      f[ i ] -= ( f[ i + 1 ] * q[k] ) ;
    }
  }
  return;
}
