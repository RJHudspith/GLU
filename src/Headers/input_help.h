/*
Copyright 2013-2025 Renwick James Hudspith

    This file (input_help.h) is part of GLU.

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
   @file input_help.h
   @brief function definition describing the possible options available to the end user
 */
#ifndef GLU_INPUT_HELP_H
#define GLU_INPUT_HELP_H

/**
   @fn void GLU_helps_those_who_help_themselves( const char *help_str )
   @brief provided an argument --help={input file param} gives available options
   @param help_str :: a single string of the form --help={something} or --autoin={something}
   @return #GLU_SUCCESS
 */
int
GLU_helps_those_who_help_themselves( const char *help_str ) ;

/**
   @fn int GLUsage( void )
   @brief print out to stdout how to use the code
 */
int
GLUsage( void ) ;

#endif
