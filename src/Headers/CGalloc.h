/*
    Copyright 2018 Renwick James Hudspith

    This file (CGalloc.h) is part of GLU.

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
   @fn CGalloc.h
   @brief prototype declarations for the CG allocations
 */
#ifndef GLU_CGALLOC_H
#define GLU_CGALLOC_H

/**
   @fn int allocate_temp_cg( struct gauges *G , struct CGtemps *CG , const GLU_bool FACG )
   @brief allocate temporaries for Coulomb gauge fixing
   @param gauges :: gauge transformation matrices
   @param CG :: CG temporaries
   @param FACG :: are we using the FACG algorithm?
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
allocate_temp_cg( struct gauges *G ,
		  struct CGtemps *CG ,
		  const GLU_bool FACG ) ;

/**
   @fn void free_temp_cg( struct gauges G , struct CGtemps CG , const GLU_bool FACG )
   @brief free the temporaries for Coulomb gauge fixing set by allocate_temp_cg()
   @param G :: gauge transformation matrices
   @param CG :: CG temporaries
   @param FACG :: did we use the FACG algorithm?
 */
void
free_temp_cg( struct gauges G ,
	      struct CGtemps CG ,
	      const GLU_bool FACG ) ;

/**
   @fn int allocate_temp_lg( struct CGtemps *CG , const GLU_bool FACG )
   @brief allocate the temporaries for Landau gauge fixing
   @param CG :: Conjugate gradient temporaries Landau gauge
   @param FACG :: did we use the FACG algorithm?
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
allocate_temp_lg( struct CGtemps *CG ,
		  const GLU_bool FACG ) ;

/**
   @fn void free_temp_lg( struct CGtemps CG , const GLU_bool FACG )
   @brief free the temporaries allocated by allocate_temp_lg()
   @param CG :: allocated conjugate gradient temporaries for Landau gauge
   @param FACG :: did we use the FACG algorithm?
 */
void
free_temp_lg( struct CGtemps CG ,
	      const GLU_bool FACG ) ;

#endif
