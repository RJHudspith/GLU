/*
    Copyright 2013-2016 Renwick James Hudspith

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

void
create_plans_DFT( struct fftw_stuff *FFTW ,
		  const size_t dims[ ND ] ,
		  const int ARR_SIZE ,
		  const int DIR ) ;

void
small_create_plans_DFT( struct fftw_small_stuff *FFTW ,
			const size_t dims[ ND ] ,
			const int DIR ) ;

void
clean_up_fftw( struct fftw_stuff FFTW ,
	       const int ARR_SIZE ) ;

void
small_clean_up_fftw( struct fftw_small_stuff FFTW ) ;

#endif // HAVE_FFTW_H
#endif
