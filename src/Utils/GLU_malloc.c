/*
    Copyright 2013-2016 Renwick James Hudspith

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

#if (defined HAVE_IMMINTRIN_H)
#include <immintrin.h>
#endif

#include "init.h" // init_navig()

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

// allocate the gauge fields
struct site *
allocate_lat( void )
{
  struct site *lat = NULL ;
  size_t i , flag = 0 ;
  if( GLU_malloc( (void**)&lat , ALIGNMENT , LVOLUME * sizeof( struct site ) ) != 0 ) {
    return NULL ;
  }
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    if( GLU_malloc( (void**)&lat[i].O , ALIGNMENT , ND * sizeof( GLU_complex* ) ) != 0 ) {
      flag = 1 ;
    }
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      if( GLU_malloc( (void**)&lat[i].O[mu] , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) != 0 ) {
	flag = 1 ;
      }
    }
  }
  if( flag == 1 ) return NULL ;
  init_navig( lat ) ;
  return lat ;
}

// free the lattice gauge field
void
free_lat( struct site *lat )
{
  size_t i , mu ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      if( lat[i].O[mu] != NULL ) {
	free( lat[i].O[mu] ) ;
      }
    }
    if( lat[i].O != NULL ) {
      free( lat[i].O ) ;
    }
  }
  if( lat != NULL ) {
    free( lat ) ;
  }
  return ;
}

// allocate an spt_size
struct s_site *
allocate_s_site( const size_t LENGTH1 ,
		 const size_t LENGTH2 ,
		 const size_t LENGTH3 )
{
  size_t i , flag = 0 ;
  struct s_site *lat = NULL ;
  if( GLU_malloc( (void**)&lat , ALIGNMENT , LENGTH1 * sizeof( struct s_site ) ) != 0 ) {
    return NULL ;
  }
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LENGTH1 ; i++ ) {
    if( GLU_malloc( (void**)&lat[i].O , ALIGNMENT , LENGTH2 * sizeof( GLU_complex* ) ) != 0 ) {
      flag = 1 ;
    }
    //
    size_t mu ;
    for( mu = 0 ; mu < LENGTH2 ; mu++ ) {
      if( GLU_malloc( (void**)&lat[i].O[mu] , ALIGNMENT , LENGTH3 * sizeof( GLU_complex ) ) != 0 ) {
	flag = 1 ;
      }
    }
  }
  if( flag == 1 ) return NULL ;
  return lat ;
}

// free 
void
free_s_site( struct s_site *lat ,
	     const size_t LENGTH1 ,
	     const size_t LENGTH2 ,
	     const size_t LENGTH3 )
{

  size_t i , mu ;
  for( i = 0 ; i < LENGTH1 ; i++ ) {
    for( mu = 0 ; mu < LENGTH2 ; mu++ ) {
      if( lat[i].O[mu] != NULL ) {
	free( lat[i].O[mu] ) ;
      }
    }
    if( lat[i].O != NULL ) {
      free( lat[i].O ) ;
    }
  }
  if( lat != NULL ) {
    free( lat ) ;
  }
  return ;
}
