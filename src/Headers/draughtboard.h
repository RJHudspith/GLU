/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (draughtboard.h) is part of GLU.

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
   @file draughtboard.h
   @brief prototype functions for setting the draughtboard
 */
#ifndef DRAUGHTBOARD_H
#define DRAUGHTBOARD_H

/**
   @fn void free_cb( struct draughtboard *db ) 
   @brief free the draughtboard
   @param db :: draughtboard structure
 */
void
free_cb( struct draughtboard *db ) ;

/**
   @fn void init_cb( struct draughtboard *db , const size_t LENGTH , const size_t DIR ) 
   @brief initialise the draughtboard
   @param db :: draughtboard structure
   @param LENGTH :: total length of the thing we are draughtboarding
   @param DIR :: number of directions in the geometry we are using
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
init_cb( struct draughtboard *db ,
	 const size_t LENGTH ,
	 const size_t DIR ) ;

int
test_db( struct site *lat ,
	 const struct draughtboard db ) ;

#endif
