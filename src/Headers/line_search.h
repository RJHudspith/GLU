/*
    Copyright 2018 Renwick James Hudspith

    This file (line_search.h) is part of GLU.

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
   @file line_search.h
   @brief line search routines for the Landau and Coulomb gauge fixing routines
 */
#ifndef GLU_LINE_SEARCH_H
#define GLU_LINE_SEARCH_H

/**
   @fn double approx_minimum( const size_t nmeas , const double alphas[ nmeas ] , const double functional[ nmeas ] )
   @brief finds the approximate minimum of alphas using GLU-bic splines
   @param nmeas :: number of measurements made
   @param alphas :: alphas tested
   @param functional :: the gauge functional at each alpha

   @return the alpha that approximately minimises the functional, or 0 
 */
double
approx_minimum( const size_t nmeas , 
		const double alphas[ nmeas ] ,
		const double functional[ nmeas ] ) ;

/**
   @fn void egauge_Landau( GLU_complex **gauge , const GLU_complex **in , const double alpha )
   @brief exponentiate the derivative into the matrix gauge
   @param gauge :: gauge rotation matrices for this iteration
   @param in :: derivative dA
   @param alpha :: step length for the SD/CG routines
 */
void
egauge_Landau( GLU_complex **gauge , 
	       const GLU_complex **in ,
	       const double alpha ) ;

/**
   @fn void exponentiate_gauge_CG( GLU_complex **gauge , const GLU_complex **in , const double alpha )
   @brief exponentiate the derivative into the matrix gauge for the Coulomb gauge fixing routines
   @param gauge :: gauge rotation matrices for this iteration
   @param in :: derivative dA
   @param alpha :: step length for the SD/CG routines
 */
void
exponentiate_gauge_CG( GLU_complex **gauge , 
		       const GLU_complex **in ,
		       const double alpha ) ;


/**
   @fn void line_search_Coulomb( double *red , GLU_complex **gauge , const struct draughtboard db , const struct site *lat , const GLU_complex **in , const size_t t )
   @brief perform the line search for the Coulomb gauge fixing routines
   @param red :: reduction array
   @param gauge :: accumulated gauge transform matrices
   @param db :: draughtboard
   @param lat :: lattice links
   @param in :: derivative of the gauge fields
   @param t :: timeslice index
 */
void
line_search_Coulomb( double *red ,
		     GLU_complex **gauge ,
		     const struct draughtboard db ,
		     const struct site *lat ,
		     const GLU_complex **in ,
		     const size_t t ) ;

/**
   @fn void line_search_Landau( double *red , GLU_complex **gauge , const struct site *lat , const GLU_complex **in )
   @brief perform the line search for Landau gauge fixing
   @param red :: reduction array
   @param gauge :: this iteration's gauge transformations
   @param lat :: lattice links
   @param in :: derivative of the gauge fields
 */
void
line_search_Landau( double *red ,
		    GLU_complex **gauge , 
		    const struct site *lat ,
		    const GLU_complex **in ) ;

#endif
