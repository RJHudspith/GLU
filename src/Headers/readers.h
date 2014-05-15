/*
    Copyright 2013 Renwick James Hudspith

    This file (readers.h) is part of GLU.

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
   @file readers.h
   @brief configuration file readers
 */
#ifndef GLU_READERS_H
#define GLU_READERS_H

/**
   @fn uint32_t lattice_reader_suNC( struct site *__restrict lat , FILE *__restrict in , const struct head_data HEAD_DATA )
   @brief reads in a NERSC configuration ...
   @param lat :: lattice gauge field
   @param in :: NERSC configuration being read
   @param HEAD_DATA :: the header data

   @returns #GLU_SUCCESS or #GLU_FAILURE
 **/
uint32_t
lattice_reader_suNC( struct site *__restrict lat , 
		     FILE *__restrict in , 
		     const struct head_data HEAD_DATA ) ;

/**
   @fn uint32_t lattice_reader_suNC_cheaper( struct site *__restrict lat , FILE *__restrict in , const struct head_data HEAD_DATA )
   @brief reads in a NERSC configuration, cheaper and a little slower than the other needed for the Q.
   @param lat :: lattice gauge field
   @param in :: NERSC configuration being read
   @param HEAD_DATA :: the header data

   @returns #GLU_SUCCESS or #GLU_FAILURE
 **/
uint32_t
lattice_reader_suNC_cheaper( struct site *__restrict lat , 
			     FILE *__restrict in , 
			     const struct head_data HEAD_DATA ) ;

#endif
