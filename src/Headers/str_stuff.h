/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (str_stuff.h) is part of GLU.

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
   @file str_stuff.h
   @brief prototype declarations for useful string functions
 */
#ifndef GLU_STR_STUFF_H
#define GLU_STR_STUFF_H

/**
   @fn void append_char( char **str , const char *tmp )
   @brief appends tmp to str, reallocating str
 */
void
append_char( char **str ,
	     const char *tmp ) ;

/**
   @fn int are_equal( const char *str_1 , const char *str_2 )
   @brief checks if two strings are the same, returns true if they are
   @return 1 or !1 if the are the same or not
 */
int
are_equal( const char *str_1 ,
	   const char *str_2 ) ;

#endif
