/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (smear.h) is part of GLU.

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
   @file smear.h
   @brief function prototypes for simple smearing algorithms here, calls staples in staples.h
 */
#ifndef GLU_SMEAR_H
#define GLU_SMEAR_H

/**
   @fn int smear3D( struct site *__restrict lat , const size_t smiters , const smearing_types type )
   @brief ND-1 link smearing without blocking transforms
   @param lat :: lattice gauge fields 
   @param smiters :: number of smearing iterations
   @param type :: type of smearing projection

   calls staples() and projectors.h. Heavily.
   @return #GLU_SUCCESS or #GLU_FAILURE 
 **/
int
smear3D( struct site *__restrict lat , 
	 const size_t smiters , 
	 const smearing_types type ) ;

/**
   @fn int smear4D( struct site *__restrict lat , const size_t smiters , const smearing_types type )
   @brief ND link smearing without blocking transforms
   @param lat :: lattice gauge fields 
   @param smiters :: number of smearing iterations
   @param type :: type of smearing projection

   calls staples() and projectors.h . Heavily.
   @return #GLU_SUCCESS or #GLU_FAILURE 
 **/
int 
smear4D( struct site *__restrict lat ,
	 const size_t smiters , 
	 const smearing_types type ) ;

#endif
