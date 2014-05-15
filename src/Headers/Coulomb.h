/*
    Copyright 2013 Renwick James Hudspith

    This file (Coulomb.h) is part of GLU.

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
   @file Coulomb.h
   @brief Coulomb gauge fixing code

   Uses either the fourier-accelerated code or the normal steepest descent.
   Performs its iterations on a slice-by-slice basis.
 */

#ifndef GLU_COULOMB_H
#define GLU_COULOMB_H

/**
   @fn int Coulomb( struct site *__restrict lat , const double accuracy , const int iter ) 
   @brief Coulomb gauge fixing
   
   @param lat :: Lattice fields
   @param accuracy :: Gauge fixing accuracy we are iterating to
   @param iter :: Maximum number of iterations before restarting

   @return #GLU_FAILURE or #GLU_SUCCESS
 */
int 
Coulomb( struct site *__restrict lat , 
	 const double accuracy , 
	 const int iter ) ;

#endif
