/*
    Copyright 2013 Renwick James Hudspith

    This file (log_derivs.h) is part of GLU.

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
   @file log_derivs.h
   @brief function definitions for the lattice derivatives used
 */

#ifndef GLU_LOG_DERIVS_H
#define GLU_LOG_DERIVS_H

/**
   @fn double log_deriv( GLU_complex sum[ HERMSIZE ] , const struct site *__restrict lat , const size_t i , const size_t MAX_DIR )
   @brief The logarithmic lattice derivative
   @param sum :: The sum of the derivative \f$ \partial_\mu A_\mu(x) \f$
   @param functional :: evaluation of the log-functional
   @param lat :: The lattice field \f$ U_\mu(x) = e^{iaA_\mu(x)} \f$
   @param i :: The site index.
   @param MAX_DIR :: The number of directions we take the derivative of. 

   @return returns the accuracy which is defined as <br>
   \f$  Tr\left( | \partial_\mu A_\mu(x) |^2 \right) \f$
 **/
double
log_deriv( GLU_complex sum[ HERMSIZE ] , 
	   double *functional ,
	   const struct site *__restrict lat , 
	   const size_t i , 
	   const size_t MAX_DIR ) ;

/**
   @fn double log_deriv_nn( GLU_complex sum[ HERMSIZE ] , const struct site *__restrict lat , const size_t i , const size_t MAX_DIR )
   @brief The (nearest neighbor) logarithmic lattice derivative
   @param sum :: The sum of the derivative \f$ \partial_\mu A_\mu(x) \f$
   @param lat :: The lattice field \f$ U_\mu(x) = e^{iaA_\mu(x)} \f$
   @param i :: The site index.
   @param MAX_DIR :: The number of directions we take the derivative of. 

   @return returns the accuracy which is defined as <br>
   \f$  Tr\left( | \partial_\mu A_\mu(x) |^2 \right) \f$
 **/
double
log_deriv_nn( GLU_complex sum[ HERMSIZE ] , 
	      const struct site *__restrict lat , 
	      const size_t i , 
	      const size_t MAX_DIR ) ;

/**
   @fn double log_deriv_nnn( GLU_complex sum[ HERMSIZE ] , const struct site *__restrict lat , const size_t i , const size_t MAX_DIR )
   @brief The (nearest neighbor) logarithmic lattice derivative
   @param sum :: The sum of the derivative \f$ \partial_\mu A_\mu(x) \f$
   @param lat :: The lattice field \f$ U_\mu(x) = e^{iaA_\mu(x)} \f$
   @param i :: The site index.
   @param MAX_DIR :: The number of directions we take the derivative of. 

   @return returns the accuracy which is defined as <br>
   \f$  Tr\left( | \partial_\mu A_\mu(x) |^2 \right) \f$
 **/
double
log_deriv_nnn( GLU_complex sum[ HERMSIZE ] , 
	       const struct site *__restrict lat , 
	       const size_t i , 
	       const size_t MAX_DIR ) ;

#endif
