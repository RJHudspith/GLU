/*
    Copyright 2013 Renwick James Hudspith

    This file (derivs.h) is part of GLU.

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
   @file lin_derivs.h
   @brief function definitions for the lattice derivatives used
 */

#ifndef GLU_LIN_DERIVS_H
#define GLU_LIN_DERIVS_H

/**
 @fn inline void constant_mul_deriv( GLU_complex sum[ HERMSIZE ] , const GLU_real constant , const GLU_complex shiftA[ HERMSIZE ] , const GLU_complex A[ HERMSIZE ] )
 @brief multiplies the derivative by a constant constant.dA
 @param sum :: derivative
 @param constant :: the constant
 @param shiftA :: A( x - \mu/2 )
 @param A :: A( x + \mu/2 )
 */
inline void
constant_mul_deriv( GLU_complex sum[ HERMSIZE ] ,
		    const GLU_real constant ,
		    const GLU_complex shiftA[ HERMSIZE ] ,
		    const GLU_complex A[ HERMSIZE ] ) ;

/**
   @fn double trace_deriv( GLU_complex *__restrict sum )
   @brief computes the trace of the derivative dA
   @param sum :: the derivative matrix
 */
double
trace_deriv( GLU_complex *__restrict sum ) ;

/**
   @fn double latt_deriv_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , const struct site *__restrict lat , const int i , const int MAX_DIR )
   @brief The lattice derivative
   @param sum :: The sum of the derivative \f$ \partial_\mu A_\mu(x) \f$
   @param lat :: The lattice field \f$ U_\mu(x) = e^{iaA_\mu(x)} \f$
   @param i :: The site index.
   @param MAX_DIR :: The number of directions we take the derivative of. 

   @return returns the accuracy which is defined as <br>
   \f$  Tr\left( | \partial_\mu A_\mu(x) |^2 \right) \f$
 **/
double
latt_deriv_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , 
			       const struct site *__restrict lat , 
			       const int i , 
			       const int MAX_DIR ) ;

/**
   @fn double fast_deriv_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , const struct site *__restrict lat , const int i )
   @brief Slightly faster version of the 4D lattice lie field derivative
   @param sum :: The sum of the derivative \f$ \partial_\mu A_\mu(x) \f$
   @param lat :: The lattice field \f$ U_\mu(x) = e^{iaA_\mu(x)} \f$
   @param sum :: Sum of the derivative of the lie field
   @param lat :: Lattice field
   @param i :: Site index

   @return returns the accuracy which is defined as <br>
   \f$  Tr\left( | \partial_\mu A_\mu(x) |^2 \right) \f$

 **/
double
fast_deriv_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , 
			       const struct site *__restrict lat , 
			       const int i ) ;

/**
   @fn double latt_derivnn_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , const struct site *__restrict lat , const int i , const int MAX_DIR )
   @brief Derivative using the stencil
   @param sum :: The sum of the derivative \f$ \partial_\mu A_\mu(x) \f$
   @param lat :: The lattice field \f$ U_\mu(x) = e^{iaA_\mu(x)} \f$
   @param i :: The site index.
   @param MAX_DIR :: The number of directions we take the derivative of. 

   @return returns the accuracy which is defined as <br>
   \f$  Tr\left( | \partial_\mu A_\mu(x) |^2 \right) \f$
 **/ 
double
latt_derivnn_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , 
				 const struct site *__restrict lat , 
				 const int i , 
				 const int MAX_DIR ) ;

/**
   @fn double fast_derivnn_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , const struct site *__restrict lat , const int i )
   @brief Derivative using the stencil, speed up significantly
   @param sum :: The sum of the derivative \f$ \partial_\mu A_\mu(x) \f$
   @param lat :: The lattice field \f$ U_\mu(x) = e^{iaA_\mu(x)} \f$
   @param i :: The site index.

   @return returns the accuracy which is defined as <br>
   \f$  Tr\left( | \partial_\mu A_\mu(x) |^2 \right) \f$
 **/ 
double
fast_derivnn_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , 
				 const struct site *__restrict lat , 
				 const int i ) ;

#endif
