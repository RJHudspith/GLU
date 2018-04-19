/*
    Copyright 2013-2018 Renwick James Hudspith

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

/**
   @fn struct site * allocate_lat( void )
   @brief allocate the lattice gauge fields
   @return the lattice gauge fields
 */
struct site *
allocate_lat( void ) ;

/**
   @fn struct s_site * allocate_spt_site( const size_t LENGTH1 , const size_t LENGTH2 , const size_t LENGTH3 )
   @brief allocate an spt_site struct lat[LENGTH1].O[LENGTH2][LENGTH3]
   @param lat :: the lattice structure
   @param LENGTH1 :: the first length
   @param LENGTH2 :: the second
   @param LENGTH3 :: the third
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
struct s_site *
allocate_s_site( const size_t LENGTH1 ,
		 const size_t LENGTH2 ,
		 const size_t LENGTH3 ) ;

/**
   @fn void free_spt_site( struct s_site *lat , const size_t LENGTH1 , const size_t LENGTH2 , const size_t LENGTH3 )
   @brief free the spt_size
*/
void
free_s_site( struct s_site *lat ,
	     const size_t LENGTH1 ,
	     const size_t LENGTH2 ,
	     const size_t LENGTH3 ) ;

/**
   @fn void free_lat( struct site *lat )
   @brief free the lattice gauge fields
   @param lat :: the lattice gauge fields
 */
void
free_lat( struct site *lat ) ;

#endif
