/*
    Copyright 2013 Renwick James Hudspith

    This file (HYP_4D.h) is part of GLU.

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
   @file HYP_4D.h
   @brief function definition for the medium memory-intensiveness code
 */
#ifndef GLU_HYP_4D_H
#define GLU_HYP_4D_H

/**
   @fn int HYPSLsmear4D( struct site *__restrict lat , const size_t smiters , const int type )
   @brief Performs a Hypecubically-blocked smearing transformation
   @param lat :: lattice gauge field
   @param smiters :: number of smearing iterations
   @param type :: type of smearing operation ( SM_LOG , SM_APE , SM_STOUT )
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
HYPSLsmear4D( struct site *__restrict lat , 
	      const size_t smiters , 
	      const int type ) ;

#endif
