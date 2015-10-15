/*
    Copyright 2013 Renwick James Hudspith

    This file (GLU_malloc.c) is part of GLU.

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
   @file corr_malloc.c
   @brief memory allocation wrapper
 */
#include "Mainfile.h"

// memalign wrapper
int 
GLU_malloc( void **memptr , 
	    const size_t alignment , 
	    const size_t size )
{
#if (defined HAVE_IMMINTRIN_H)
  return posix_memalign( memptr , alignment , size ) ;
#else
  *memptr = malloc( size ) ;
  return ( *memptr == NULL ) ? 1 : 0 ;
#endif
}
