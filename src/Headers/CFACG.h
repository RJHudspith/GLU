/*
Copyright 2013-2025 Renwick James Hudspith

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
   @fn size_t Coulomb_FA( struct site  *__restrict lat , struct fftw_stuff *FFTW , const double accuracy , const size_t max_iter , const GLU_bool FACG ) ;
   @brief Coulomb gauge fixing codes
   @param lat :: lattice links
   @param FFTW :: fftw temporaries and plans
   @param accuracy :: the accuracy we wish to attain
   @param max_iter :: the maximum number of (per slice) iterations we wish to have before random transform
   @param FACG :: are we using the CG routines?
   @return the maximum number of iterations (i.e. the sum of each slice's) or 0 if something went wrong
 */
size_t
Coulomb_FA( struct site  *__restrict lat , 
	    struct fftw_stuff *FFTW ,
	    const double accuracy ,
	    const size_t max_iter ,
	    const GLU_bool FACG ) ;

#endif
