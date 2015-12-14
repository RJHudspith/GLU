/*
    Copyright 2013 Renwick James Hudspith

    This file (CFACG.h) is part of GLU.

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
   @file CFACG.h
   @brief Fourier Accelerated Coulomb gauge fixing routines
 */
#ifndef GLU_CFACG_H
#define GLU_CFACG_H

/**
   @fn size_t Coulomb_FACG( struct site  *__restrict lat , GLU_complex  *__restrict *__restrict out , GLU_complex  *__restrict *__restrict in , const void *__restrict forward , const void *__restrict backward , const GLU_real * __restrict p_sq , const double accuracy , const size_t max_iter ) ;
   @brief Conjugate gradient Fourier Accelerated Coulomb gauge fixing code
   @param lat :: lattice links
   @param out :: FFTW temporary
   @param in :: derivative of the fields
   @param forward :: FFTW forward transform plan
   @param backward :: FFTW backward transform plan
   @param accuracy :: the accuracy we wish to attain
   @param max_iter :: the maximum number of (per slice) iterations we wish to have before random transform
   @return the maximum number of iterations (i.e. the sum of each slice's) or 0 if something went wrong
 */
size_t
Coulomb_FACG( struct site  *__restrict lat , 
	      GLU_complex  *__restrict *__restrict out , 
	      GLU_complex  *__restrict *__restrict in , 
	      const void *__restrict forward , 
	      const void *__restrict backward , 
	      const GLU_real * __restrict p_sq ,
	      const double accuracy ,
	      const size_t max_iter ) ;
/**
   @fn size_t Coulomb_FASD( struct site  *__restrict lat , GLU_complex  *__restrict *__restrict out , GLU_complex  *__restrict *__restrict in , const void *__restrict forward , const void *__restrict backward , const GLU_real * __restrict p_sq , const double accuracy , const size_t max_iter ) ;
   @brief steepest descent Fourier Accelerated Coulomb gauge fixing code
   @param lat :: lattice links
   @param out :: FFTW temporary
   @param in :: derivative of the fields
   @param forward :: FFTW forward transform plan
   @param backward :: FFTW backward transform plan
   @param psq :: precomputed psq table
   @param accuracy :: the accuracy we wish to attain
   @param max_iter :: the maximum number of (per slice) iterations we wish to have before random transform
   @return the maximum number of iterations (i.e. the sum of each slice's) or 0 if something went wrong
 */
size_t
Coulomb_FASD( struct site  *__restrict lat , 
	      GLU_complex  *__restrict *__restrict out , 
	      GLU_complex  *__restrict *__restrict in , 
	      const void *__restrict forward , 
	      const void *__restrict backward , 
	      const GLU_real * __restrict p_sq ,
	      const double accuracy ,
	      const size_t max_iter ) ;

/**
   @fn void query_probes_Coulomb( void )
   @brief prints to stdout the probe alphas we use
 */
void
query_probes_Coulomb( void ) ;

#endif
