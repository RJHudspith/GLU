/*
    Copyright 2013 Renwick James Hudspith

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
   @fn void query_probes_Landau( void ) 
   @brief tells us the probe alphas we have set
   prints to stdout
 */
void
query_probes_Landau( void ) ;

/**
   @fn size_t FACG( struct site *__restrict lat , GLU_complex *__restrict *__restrict gauge , GLU_complex *__restrict *__restrict out , GLU_complex *__restrict *__restrict in , const void *__restrict forward , const void *__restrict backward , const GLU_real *__restrict p_sq , double *th , const double acc , const size_t max_iters )
   @brief Fourier Accelerated CG Landau gauge fixing
   @param lat :: the gauge links (overwritten)
   @param gauge :: gauge transformation matrices (overwritten)
   @param out :: FFT'd array of size [HERMSIZE][LVOLUME]
   @param in :: FFT'd array of size [HERMSIZE][LVOLUME]
   @param forward :: forward FFTW plan
   @param backward :: backward FFTW plan
   @param psq :: psq array, rescaled by 1/V
   @param th :: theta, the quantity that controls when we stop
   @param acc :: the targeted gauge fixing accuracy
   @param max_iters :: the maximum number of iterations we will perform
   @warning is in general, faster than the steepest descent method at the cost of a little extra storage
 */
size_t
FACG( struct site *__restrict lat ,
      GLU_complex *__restrict *__restrict gauge ,
      GLU_complex *__restrict *__restrict out ,
      GLU_complex *__restrict *__restrict in ,
      const void *__restrict forward ,
      const void *__restrict backward ,
      const GLU_real *__restrict p_sq ,
      double *th ,
      const double acc ,
      const size_t max_iters ) ;

/**
   @fn int FASD( struct site *__restrict lat , GLU_complex *__restrict *__restrict gauge , GLU_complex *__restrict *__restrict out , GLU_complex *__restrict *__restrict in , const void *__restrict forward , const void *__restrict backward , const GLU_real *__restrict p_sq , double *th , const double acc , const size_t max_iters )
   @brief Fourier Accelerated CG Landau gauge fixing
   @param lat :: the gauge links (overwritten)
   @param gauge :: gauge transformation matrices (overwritten)
   @param out :: FFT'd array of size [HERMSIZE][LVOLUME]
   @param in :: FFT'd array of size [HERMSIZE][LVOLUME]
   @param forward :: forward FFTW plan
   @param backward :: backward FFTW plan
   @param psq :: psq array, rescaled by 1/V
   @param th :: theta, the quantity that controls when we stop
   @param acc :: the targeted gauge fixing accuracy
   @param max_iters :: the maximum number of iterations we will perform
   @warning probably deprecated for the CG
 */
size_t
FASD( struct site *__restrict lat ,
      GLU_complex *__restrict *__restrict gauge ,
      GLU_complex *__restrict *__restrict out ,
      GLU_complex *__restrict *__restrict in ,
      const void *__restrict forward ,
      const void *__restrict backward ,
      const GLU_real *__restrict p_sq ,
      double *th ,
      const double acc ,
      const size_t max_iters ) ;


/**
   @fn size_t FASD_SMEAR( struct site *__restrict lat , GLU_complex *__restrict *__restrict gauge , GLU_complex *__restrict *__restrict out , GLU_complex *__restrict *__restrict in , const void *__restrict forward , const void *__restrict backward , const GLU_real *__restrict p_sq , double *th , const double acc , const size_t max_iters )
   @brief Fourier Accelerated CG Landau gauge fixing
   @param lat :: the gauge links (overwritten)
   @param gauge :: gauge transformation matrices (overwritten by the product at each step)
   @param out :: FFT'd array of size [HERMSIZE][LVOLUME]
   @param in :: FFT'd array of size [HERMSIZE][LVOLUME]
   @param forward :: forward FFTW plan
   @param backward :: backward FFTW plan
   @param psq :: psq array, rescaled by 1/V
   @param th :: theta, the quantity that controls when we stop
   @param acc :: the targeted gauge fixing accuracy
   @param max_iters :: the maximum number of iterations we will perform
 */
size_t 
FASD_SMEAR( struct site *__restrict lat ,
	    GLU_complex *__restrict *__restrict gauge ,
	    GLU_complex *__restrict *__restrict out ,
	    GLU_complex *__restrict *__restrict in ,
	    const void *__restrict forward ,
	    const void *__restrict backward ,
	    const GLU_real *__restrict p_sq ,
	    double *th ,
	    const double acc ,
	    const size_t max_iters ) ;

#endif
