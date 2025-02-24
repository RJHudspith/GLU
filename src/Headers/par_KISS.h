/*
Copyright 2013-2025 Renwick James Hudspith

    This file (par_KISS.h) is part of GLU.

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
   @file par_KISS.h
   @brief parallel KISS rng
 */
#ifndef GLU_PAR_KISS_H
#define GLU_PAR_KISS_H

/**
  @fn void free_par_KISS( void )
  @brief free the KISS rng
 */
void
free_par_KISS( void ) ;

/**
   @fn par_KISS_dbl( const uint32_t thread )
   @brief redturn a uniformly distributed double value
 */
double
par_KISS_dbl( const uint32_t thread ) ;

/**
   @fn int read_par_KISS_table( FILE *rng_file ) 
   @brief read the KISS table from a file
   @param rng_file :: the file with the rng state in
 */
int
read_par_KISS_table( FILE *rng_file ) ;

/**
   @fn int write_par_KISS_table( FILE *rng_file )
   @brief write the KISS rng state to a file
   @param rng_file :: the file being written to
 */
int
write_par_KISS_table( FILE *rng_file ) ;

#endif
