/*
    Copyright 2013 Renwick James Hudspith

    This file (HIREP.h) is part of GLU.

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
   @file HIREP.h
   @brief function defs for IO with HiRep configuration files
 */

#ifndef GLU_HIREP_H
#define GLU_HIREP_H

/**
   @fn int read_gauge_field( struct site *__restrict lat , FILE *__restrict in , uint32_t *chksum )
   @brief reads in a configuration from a HiRep format file 
   @param lat :: lattice gauge fields
   @param in :: file to be read in
   @param chksum :: compute the checksum of the data for no reason whatsoever
   @return #GLU_SUCCESS or #GLU_FAILURE
**/
int
read_gauge_field( struct site *__restrict lat , 
		  FILE *__restrict in , 
		  uint32_t *chksum ) ;

/**
   @fn void write_gauge_field( const struct site *__restrict lat , FILE *__restrict outfile )
   @brief writes out a configuration file in the HiRep order
   @param lat :: lattice gauge field
   @param outfile :: file we are outputting to
 **/
void
write_gauge_field( const struct site *__restrict lat ,
		   FILE *__restrict outfile ) ;

#endif
