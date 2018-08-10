/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (MAG.h) is part of GLU.

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
   @file MAG.h
   @brief fixing to the Maximal Axial Gauge, Axial gauge, or residual (after Coulomb) is performed here
 */
#ifndef GLU_MAG_H
#define GLU_MAG_H

/**
   @fn void axial_gauge( struct site *lat , const size_t DIR )
   @brief Fixes to the axial gauge i.e \f$ A_{DIR} = 0. \f$
   @param lat :: lattice fields
   @param DIR :: direction to set fields to unity

   \f[
   A_{DIR}(x) = 0.
   \f]
   @return #GLU_SUCCESS or #GLU_FAILURE
 **/
int
axial_gauge( struct site *lat ,
	     const size_t DIR ) ;

/**
   @fn void mag( struct site *lat )
   @brief wrapper for the maximal axial gauge (MAG) fixing
   @param lat :: lattice fields   
   @warning has a reunitarization method included, to sharpen the links a little.
   @return #GLU_SUCCESS or #GLU_FAILURE
 **/
int 
mag( struct site *lat ) ;


/**
   @fn void residual_fix( struct site *lat )
   @brief fixes the extra degree of freedom left from fixing to Coulomb gauge
   @param lat :: the lattice links
   This function computes
   \f[
   U_t(x,t) = g(t) U_t(x,t) g^{\dagger}(t+1) \qquad U_i(x,t) = g(t) U_i(x,t) g^{\dagger}(t)
   \f]
   <\br>
   where the gauge transformation matrices are defined as
   \f[
   g(t) = g(t-1)\frac{1}{\prod_{i<ND-1}L_i}\textrm{Proj}_{SU(NC)}\left( \sum_x U_t(x,t) \right) \qquad g(0) = I
   \f]

   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
residual_fix( struct site *lat ) ;

#endif
