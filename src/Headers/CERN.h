/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (CERN.h) is part of GLU.

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
#ifndef GLU_CERN_H
#define GLU_CERN_H

/**
   @fn int read_CLS_field( struct site *lat , FILE *in , uint32_t *chksum )
   @param lat :: lattice gauge links
   @param in :: infile to be read
   @param chksum :: dummy parameter not used
   @brief read in a CERN gauge field
 */
int
read_CLS_field( struct site *lat , 
		FILE *in , 
		uint32_t *chksum ) ;

/**
   @fn void write_CLS_field( const struct site *lat ,  FILE *outfile )
   @param lat :: lattice gauge field
   @param outfile :: string we are writing out to
   @brief write a CERN gauge field
 */
void
write_CLS_field( const struct site *lat ,
		 FILE *outfile ) ;

#endif
