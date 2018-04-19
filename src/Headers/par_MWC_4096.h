/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (par_MWC_4096.h) is part of GLU.

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
   @file par_MWC_4096.h
   @brief prototype declarations for the MWC 4096
 */
#ifndef GLU_PAR_MWC_4096_H
#define GLU_PAR_MWC_4096_H

/**
   @fn void free_par_MWC_4096( void ) 
   @brief free the RNG table
 */
void
free_par_MWC_4096( void ) ;

/**
   @fn void GLU_set_par_MWC_4096_table( const uint32_t seed[ Latt.Nthreads ] )
   @brief set the parallel table
   @param seed :: seeds for the RNG
 */
void
GLU_set_par_MWC_4096_table( const uint32_t seed[ Latt.Nthreads ] ) ;

/**
   @fn double par_MWC_4096_dbl( const uint32_t thread )
   @brief create a double precision random number
   @param thread :: parallel thread index
 */
double 
par_MWC_4096_dbl( const uint32_t thread ) ;

/**
   @fn int read_par_MWC_4096_table( FILE *rng_file )
   @brief read the Well RNG table from rng_file
   @param rng_file :: file to read the table from
   @return #GLU_SUCCSS or #GLU_FAILURE
 */
int
read_par_MWC_4096_table( FILE *rng_file ) ;

/**
   @fn void write_par_MWC_4096_table( FILE *rng_file )
   @brief write the WELL RNG table to rng_file
   @param rng_file :: file to write out the table
 */
void
write_par_MWC_4096_table( FILE *rng_file ) ;

#endif
