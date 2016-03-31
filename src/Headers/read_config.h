/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (read_config.h) is part of GLU.

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
   @file read_config.h
   @brief prototype functions for header reading and checksum calculating
 */
#ifndef GLU_READ_CONFIG
#define GLU_READ_CONFIG

/**
   @fn int checks( struct site *__restrict lat , uint32_t chksum , struct head_data HEAD_DATA )
   @brief perform a calculation/check of some checksums against the header ...
   @param lat :: lattice gauge field
   @param chksum :: the checksum
   @param HEAD_DATA :: the header data we compare to

   @warning I am optimistic here, I only claim we are a #GLU_FAILURE if all the checksums fail
   @returns the #GLU_SUCCESS or #GLU_FAILURE
 **/
int 
checks( struct site *__restrict lat , 
	uint32_t chksum ,
	struct head_data HEAD_DATA ) ;

/**
   @fn int get_config_SUNC( FILE *__restrict CONFIG , struct site *__restrict lat , const struct head_data HEAD_DATA )
   @brief wrapper for reading a configuration
   @param CONFIG :: our configuration file
   @param lat :: the lattice gauge field being read
   @param HEAD_DATA :: uses the header data
   @returns #GLU_SUCCESS or #GLU_FAILURE
 **/
int
get_config_SUNC( FILE *__restrict CONFIG , 
		 struct site *__restrict lat ,
		 const struct head_data HEAD_DATA ) ;

#endif
