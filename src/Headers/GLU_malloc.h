/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (GLU_malloc.h) is part of GLU.

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
   @file GLU_malloc.h
   @brief prototype declarations for memory allocations in GLU
 */
#ifndef GLU_MALLOC_H
#define GLU_MALLOC_H

/**
   @fn int GLU_malloc( void **memptr , const size_t alignment , const size_t size )
   @brief memory allocation wrapper
   @param memptr :: pointer to memory
   @param alignment :: memory offset
   @param size :: number of elements in array
 */
int 
GLU_malloc( void **memptr , 
	    const size_t alignment , 
	    const size_t size ) ;

#endif
