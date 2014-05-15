/*
    Copyright 2013 Renwick James Hudspith

    This file (read_headers.h) is part of GLU.

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
   @file read_headers.h
   @brief reads the header information from the file format specified
 */
#ifndef GLU_READ_HEADERS_H
#define GLU_READ_HEADERS_H

/**
   @fn int read_header( FILE *__restrict infile , struct head_data *__restrict HEAD_DATA , const GLU_bool VERB )
   @brief read the header of the specified configuration file's header
   @param infile :: configuration file
   @param HEAD_DATA :: header information taken from the file
   @param VERB :: verbose output or not?

   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
read_header( FILE *__restrict infile ,
	     struct head_data *__restrict HEAD_DATA ,
	     const GLU_bool VERB ) ;

#endif
