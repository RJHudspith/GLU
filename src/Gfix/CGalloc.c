/*
    Copyright 2018 Renwick James Hudspith

    This file (CGalloc.c) is part of GLU.

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
   @file CGalloc.c
   @brief Gauge fixing memory allocation routines
 */
#include "Mainfile.h"      // for all the definitions

#include "draughtboard.h"  // draughtboarding

// allocates the temporaries
int
allocate_temp_cg( struct gauges *G ,
		  struct CGtemps *CG ,
		  const GLU_bool FACG )
{
  GLU_bool FLAG = GLU_FALSE ;
  size_t i ;
  if( ( CG -> red = calloc( CLINE * Latt.Nthreads , sizeof( double ) ) ) == NULL ) {
    FLAG = GLU_TRUE ;  goto end ;
  }
  // checkerboard for the line search
  if( init_cb( &CG -> db , LCU , ND-1 ) == GLU_FAILURE ) {
    FLAG = GLU_TRUE ;  goto end ;
  }
  // temporary field allocations
  if( GLU_malloc( (void**)&G->g      , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&G->g_end  , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&G->g_up   , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ) {
    FLAG = GLU_TRUE ; goto end ;
  }
  // FACG has two other temporaries
  if( FACG == GLU_TRUE ) {
    if( GLU_malloc( (void**)&CG -> sn     , ALIGNMENT , TRUE_HERM * sizeof( GLU_complex* ) ) != 0 ||
	GLU_malloc( (void**)&CG -> in_old , ALIGNMENT , TRUE_HERM * sizeof( GLU_complex* ) ) != 0 ) {
      FLAG = GLU_TRUE ; goto end ;
    }
  }
   // reduce some overhead here
#pragma omp parallel
  {
    if( FACG == GLU_TRUE ) {
      #pragma omp for private(i)
      for( i = 0 ; i < TRUE_HERM ; i ++  ) {
	if( GLU_malloc( (void**)&CG -> sn[i] , ALIGNMENT , LCU * sizeof( GLU_complex ) ) != 0 ||
	    GLU_malloc( (void**)&CG -> in_old[i] , ALIGNMENT , LCU * sizeof( GLU_complex ) ) != 0 ) {
	  FLAG = GLU_TRUE ;
	}
	size_t j ;
	for( j = 0 ; j < LCU ; j++ ) {
	  CG -> sn[i][j] = CG -> in_old[i][j] = 0.0 ;
	}
      }
    }
    #pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      // allocations
      if( GLU_malloc( (void**)&G->g[i]     , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) != 0 ||
	  GLU_malloc( (void**)&G->g_up[i]  , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) != 0 ||
	  GLU_malloc( (void**)&G->g_end[i] , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) != 0 ) {
	FLAG = GLU_TRUE ;
      }
      if( FLAG != GLU_TRUE ) {
	// set gauge matrices to the identity
	identity( G->g_up[i] ) ;
	identity( G->g[i] ) ;
	identity( G->g_end[i] ) ;
      }
    }
  }
 end :
  if( FLAG == GLU_TRUE ) {
    return GLU_FAILURE ;
  } else {
    return GLU_SUCCESS ;
  }
}

// allocates the temporaries
int
allocate_temp_lg( struct CGtemps *CG ,
		  const GLU_bool FACG )
{
  CG -> sn = NULL ; CG -> in_old = NULL ; CG -> red = NULL ;
  size_t i ;
  GLU_bool FLAG = GLU_FALSE ;
  if( ( CG -> red = calloc( CLINE * Latt.Nthreads , sizeof( double ) ) ) == NULL ) {
    FLAG = GLU_TRUE ;  goto end ;
  }

  if( FACG == GLU_FALSE ) goto end ;
  
  if( GLU_malloc( (void**)&CG->sn     , ALIGNMENT , TRUE_HERM * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&CG->in_old , ALIGNMENT , TRUE_HERM * sizeof( GLU_complex* ) ) != 0 ) {
    FLAG = GLU_TRUE ;  goto end ;
  }

  // allocate the conjugate directions
  #pragma omp parallel for private(i)
  for( i = 0 ; i < TRUE_HERM ; i++ ) {
    if( GLU_malloc( (void**)&CG->sn[i]     , ALIGNMENT , LVOLUME * sizeof( GLU_complex ) ) != 0  ||
	GLU_malloc( (void**)&CG->in_old[i] , ALIGNMENT , LVOLUME * sizeof( GLU_complex ) ) != 0 ) {
      FLAG = GLU_TRUE ;
    }
  }
 end:
  if( FLAG == GLU_TRUE ) {
    return GLU_FAILURE ;
  } else {
    return GLU_SUCCESS ;
  }
}

// free Coulomb temporaries
void
free_temp_cg( struct gauges G ,
	      struct CGtemps CG ,
	      const GLU_bool FACG )
{
  size_t i ;
  // free the rotated temporary
  if( CG.red != NULL ) {
    free( CG.red ) ;
  }
  // free the CG stuff
  if( FACG == GLU_TRUE ) {
    if( CG.sn != NULL ) {
      for( i = 0 ; i < TRUE_HERM ; i++ ) {
	free( CG.sn[i] ) ;
      }
    }
    if( CG.in_old != NULL ) {
      for( i = 0 ; i < TRUE_HERM ; i++ ) {
	free( CG.in_old[i] ) ;
      }
    }
    free( CG.in_old ) ;
    free( CG.sn ) ;
  }
  // free the draughtboard
  free_cb( &CG.db ) ;
  // free temporary gauge transformation matrices
  if( G.g != NULL ) {
#pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i ++  ) {
      free( G.g[i]     ) ; 
    }
  }
  if( G.g_up != NULL ) {
#pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i ++  ) {
      free( G.g_up[i]  ) ; 
    }
  }
  if( G.g_end != NULL ) {
#pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i ++  ) {
      free( G.g_end[i] ) ; 
    }
  }
  free( G.g     ) ; 
  free( G.g_up  ) ; 
  free( G.g_end ) ;
  return ;
}

// free landau temporaries
void
free_temp_lg( struct CGtemps CG ,
	      const GLU_bool FACG )
{
  size_t i ;
  if( CG.red != NULL ) {
    free( CG.red ) ;
  }
  // free the CG stuff
  if( FACG == GLU_TRUE ) {
    if( CG.sn != NULL ) {
      for( i = 0 ; i < TRUE_HERM ; i++ ) {
	free( CG.sn[i] ) ;
      }
    }
    if( CG.in_old != NULL ) {
      for( i = 0 ; i < TRUE_HERM ; i++ ) {
	free( CG.in_old[i] ) ;
      }
    }
    free( CG.in_old ) ;
    free( CG.sn ) ;
  }
  return ;
}
