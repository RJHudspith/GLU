/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (writers.h) is part of GLU.

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
   @file writers.h
   @brief prototype function to write out lattice configuration
 */
#ifndef GLU_WRITERS_H
#define GLU_WRITERS_H

/**
   @fn int write_lat ( struct site *__restrict lat , FILE *__restrict out , const GLU_output type , char *__restrict details )
   @brief writes out our lattice data
   @param lat :: lattice gauge field
   @param out :: configuration file being written
   @param type :: output type specified, should be #config_size
   @param details :: string of information to put into the "INFO" of the NERSC header
 */
int
write_lat ( struct site *__restrict lat , 
	    FILE *__restrict out , 
	    const GLU_output type , 
	    const char *__restrict details ) ;

#endif
