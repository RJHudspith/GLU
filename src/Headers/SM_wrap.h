/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (SM_wrap.h) is part of GLU.

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
   @file SM_wrap.h
   @brief wrapper for the smearing transformations
 */
#ifndef GLU_SM_WRAP_H
#define GLU_SM_WRAP_H

/**
   @fn int SM_wrap_struct( struct site *__restrict lat , const struct sm_info SMINFO )
   @brief this is the smearing-wrapping function.
   @param lat :: lattice gauge field
   @param SMINFO :: smearing information <br>

   in here you will find <br>
   APE , HYP , STOUT , HEX , LOG , HYL and Wilson flow #smearing_types
   in ND and ND-1 dimensions <br>
   returns #GLU_SUCCESS or #GLU_FAILURE
 **/
int 
SM_wrap_struct( struct site *__restrict lat ,
		const struct sm_info SMINFO ) ;

#endif 
