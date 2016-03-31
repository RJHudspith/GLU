/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (expMat.h) is part of GLU.

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
   @file expMat.h
   @brief prototype functions for expMat.c
 */

#ifndef GLU_EXPMAT_H
#define GLU_EXPMAT_H

#if !( defined HAVE_LAPACKE_H || defined HAVE_GSL ) && ( NC > 3 )
/**
   @fn void init_factorial( void )
   @brief allocates the memory of our (inverse) factorial for our brute force matrix exponentiation routine
 */
void
init_factorial( void ) ;

/**
   @fn void free_factorial( void )
   @brief frees the precomputed factorial for the brute-force exponentiation routines
 */
void
free_factorial( void ) ;
#endif

/**
   @fn void exponentiate( GLU_complex U[ NCNC ] , const GLU_complex Q[ NCNC ] )
   @brief Our old exact exponentiation routine of Q into U. e.g. \f$ U_\mu(x) = e^{Q_\mu(x)}\f$
   @param U :: Link matrix
   @param Q :: Lie matrix
   MP exponentiation using the full (redundant) Q uses inlined eigenvalue and f calculations instead of those in solver.h and effs.h
 **/
void 
exponentiate( GLU_complex U[ NCNC ] , 
	      const GLU_complex Q[ NCNC ] ) ;

/**
   @fn void exponentiate_short( GLU_complex U[ NCNC ] , const GLU_complex Q[ HERMSIZE ] )
   @brief Our old exact exponentiation routine of Q into U. e.g. \f$ U_\mu(x) = e^{Q_\mu(x)}\f$
   @param U :: Link matrix
   @param Q :: Lie matrix
   MP exponentiation using the shortened version of Q. uses inlined eigenvalue and f calculations instead of those in solver.h and effs.h
 **/
void
exponentiate_short( GLU_complex U[ NCNC ] , 
		    const GLU_complex Q[ HERMSIZE ] ) ;

#endif
