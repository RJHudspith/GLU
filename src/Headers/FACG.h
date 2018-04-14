/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (FACG.h) is part of GLU.

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
   @file FACG.h
   @brief protype functions for the Fourier Accelerated Conjugate Gradient Landau gauge fixing routine
 */
#ifndef GLU_FACG_H
#define GLU_FACG_H

/**
 */
size_t
FACG( struct site *__restrict lat ,
      GLU_complex *__restrict *__restrict gauge ,
      struct fftw_stuff FFTW ,
      double *th ,
      const double acc ,
      const size_t max_iters ) ;

/**
 */
size_t
FASD( struct site *__restrict lat ,
      GLU_complex *__restrict *__restrict gauge ,
      struct fftw_stuff FFTW ,
      double *th ,
      const double acc ,
      const size_t max_iters ) ;


/**

 */
size_t 
FASD_SMEAR( struct site *__restrict lat ,
	    GLU_complex *__restrict *__restrict gauge ,
	    struct fftw_stuff FFTW ,
	    double *th ,
	    const double acc ,
	    const size_t max_iters ) ;

#endif
