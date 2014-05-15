/*
    Copyright 2013 Renwick James Hudspith

    This file (ND_generic_HYP.h) is part of GLU.

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
   @file ND_generic_HYP.h
   @brief prototype functions for the generic HYP code
   @ingroup Smear
 */

#ifndef GLU_ND_GENERIC_HYP_H
#define GLU_ND_GENERIC_HYP_H

/**
   @fn void HYsmearND( struct site *__restrict lat , const int smiters , const int type , const int directions )
   @brief computes the ND-generic blocking transform via slow recursion
   @param lat :: lattice gauge field
   @param smiters :: number of smearing iterations being performed
   @param type :: smearing/projection type
   @param directions :: maximum dimensions of the problem
   Recursion was unfortunately the only way I could think of doing this.
   @warning this code is incredibly slow
 */
void 
HYsmearND( struct site *__restrict lat , 
	   const int smiters , 
	   const int type ,
	   const int directions ) ;

#endif
