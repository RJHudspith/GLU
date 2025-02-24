/*
Copyright 2013-2025 Renwick James Hudspith

    This file (plan_ffts.h) is part of GLU.

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
   @file plan_ffts.h
   @brief prototype function wrappers for the ffts used 
 */
#ifndef GLU_PLAN_FFTS_H
#define GLU_PLAN_FFTS_H

#ifdef HAVE_FFTW3_H

/**
   @fn void clean_up_fftw( struct fftw_stuff FFTW , const size_t ARR_SIZE ) 
   @brief memory deallocations for fftw_stuff struct allocated in create_plans_DFT()
   @param FFTW :: collection of temporary arrays
   @param ARR_SIZE :: corresponds to ARR_SIZE in create_plans_DFT
 */
void
clean_up_fftw( struct fftw_stuff FFTW ,
	       const size_t ARR_SIZE ) ;

/**
   @fn void create_plans_DFT( struct fftw_stuff *FFTW , const size_t dims[ ND ] , const int ARR_SIZE , const size_t DIR )
   @brief create the plans and allocate arrays for the DFT
   @param FFTW :: temporary arrays
   @param dims :: dimensions of the fft
   @param ARR_SIZE :: first index length of FFTW.**in and FFTW.**out
   @param DIR :: number of dimensions of the FFT, need not be #ND
 */
void
create_plans_DFT( struct fftw_stuff *FFTW ,
		  const size_t dims[ ND ] ,
		  const size_t ARR_SIZE ,
		  const size_t DIR ) ;

/**
   @fn void small_clean_up_fftw( struct fftw_small_stuff FFTW ) 
   @brief deallocate memory and free plans created by small_create_plans_DFT()
   @param FFTW :: temporary array allocations
*/
void
small_clean_up_fftw( struct fftw_small_stuff FFTW ) ;

/**
   @fn void small_create_plans_DFT( struct fftw_small_stuff *FFTW , const size_t dims[ ND ] , const size_t DIR )
   @brief create plans and allocate arrays of FFTW.*in etc
   @param FFTW :: temporary fft arrays
   @param dims :: dimensions of the FFT in GLU order
   @param DIR :: dimensionality of the FFT, need not be #ND
 */
void
small_create_plans_DFT( struct fftw_small_stuff *FFTW ,
			const size_t dims[ ND ] ,
			const size_t DIR ) ;

#endif // HAVE_FFTW_H
#endif
