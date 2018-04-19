/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (input_reader.h) is part of GLU.

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
   @file input_reader.h
   @brief function definitions for getting information from our input file and telling the code what to do
 */
#ifndef GLU_INPUT_READER_H
#define GLU_INPUT_READER_H

/**
   @fn int get_input_data( struct infile_data *INFILE , const char *file_name ) 
   @brief read in all the information from the input file
   @param INFILE :: the struct of all the input file
   @param file_name :: the input file's name

   If any of the information isn't there or doesn't make sense the code
   will complain and probably exit.

   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
get_input_data( struct infile_data *INFILE ,
		const char *file_name ) ;

#endif
