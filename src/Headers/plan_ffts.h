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

/**
   @fn short int parallel_ffts( void )
   @brief need to initialize the openmp FFTW 
   This function calls omp_get_num_threads() and passes the
   answer to FFTW to use as the number of threads for the
   FFT.
 **/
short int
parallel_ffts( void ) ;

/**
   @fn void create_plans_DFT( fftw_plan *__restrict forward , fftw_plan *__restrict backward , GLU_complex *__restrict *__restrict in , GLU_complex *__restrict *__restrict out , const int ARR_SIZE , const int DIR )
   @brief creates an array (size ARR_SIZE) of complex to complex FFTW plans 
   @param forward :: forward FFT
   @param backward :: backward FFT
   @param in :: temporary field going in
   @param out :: temporary field out
   @param dims :: dimensions of the transform in our lexicographical ordering
   @param ARR_SIZE :: out and in array sizes
   @param DIR :: number of dimensions of the transform
   This function creates an array of plans that fourier
   transform the arrays **in and **out. ARR_SIZE is the
   size of out and in and DIR is commonly ND or ND-1.

   @warning out and in should be the same size
   <br>
   if #NOT_CONDOR_MODE is activated the functions look for wisdom in a folder called Local/Wisdom/ in the source directory
**/
void
create_plans_DFT( fftw_plan *__restrict forward , 
		  fftw_plan *__restrict backward ,
		  GLU_complex *__restrict *__restrict in , 
		  GLU_complex *__restrict *__restrict out ,
		  const size_t dims[ ND ] ,
		  const int ARR_SIZE ,
		  const int DIR ) ;

/**
   @fn void create_plans_DHT( fftw_plan *__restrict plan , GLU_real *__restrict *__restrict in , GLU_real *__restrict *__restrict out , const int ARR_SIZE , const int DIR ) 
   @brief Discrete Hartley Transform (DHT) code, used by the U(1) field generation. 
   @param plan :: the DHT plan
   @param in :: temporary field going in
   @param out :: temporary field out
   @param out :: temporary field out   
   @param dims :: dimensions of the transform
   @param ARR_SIZE :: the size of the array being FFT'd
   @param DIR :: number of dimensions of the transform
   There is no backward plan as it is its own inverse.

   @warning out and in must be the same size
   <br>
   if NOT_CONDOR_MODE is activated the functions look for wisdom in a folder called Local/Wisdom/ in the source directory
 **/
void
create_plans_DHT( fftw_plan *__restrict plan , 
		  GLU_real *__restrict *__restrict in , 
		  GLU_real *__restrict *__restrict out ,
		  const size_t dims[ ND ] ,
		  const int ARR_SIZE ,
		  const int DIR ) ;

/**
   @fn void small_create_plans_DFT( fftw_plan *__restrict forward , fftw_plan *__restrict backward , GLU_complex *__restrict in , GLU_complex *__restrict out , const int DIR )
   @brief creates a forward and a backward complex to complex fourier transform
   @param forward :: forward FFT
   @param backward :: backward FFT
   @param in :: temporary field going in
   @param out :: temporary field out
   @param dims :: dimensions of the transform
   @param DIR :: number of dimensions of the transform
   Same thing as void create_plans_DFT() but in and out are just flat arrays this is called primarily by the cutting routines and by the slow landau code in Landau.c <br>
   The forward transform is in -> out <br>
   The backward transform is out -> in <br>

   @warning out and in should be the same size
   <br>
   if #NOT_CONDOR_MODE is activated the functions look for wisdom in a folder called Local/Wisdom/ in the source directory

**/
void
small_create_plans_DFT( fftw_plan *__restrict forward , 
			fftw_plan *__restrict backward ,
			GLU_complex *__restrict in , 
			GLU_complex *__restrict out ,
			const size_t dims[ ND ] ,
			const int DIR ) ;

#endif // HAVE_FFTW_H
#endif
