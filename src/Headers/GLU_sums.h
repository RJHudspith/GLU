/*
Copyright 2013-2025 Renwick James Hudspith

    This file (GLU_sums.h) is part of GLU.

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
   @file GLU_sums.h
   @brief prototype declarations for robust summations
 */
#ifndef GLU_SUMS_H
#define GLU_SUMS_H

/**
   @fn double dandc_sum( const double *a , const size_t lo , const size_t hi )
   @brief recursive divide and conquer sum
   @param a :: array of values
   @param lo :: lower index of values (i.e. 0)
   @param hi :: length of the array @a
   @return the sum of the values in @a
   Note:: Error should go like hi.log( hi )
 */
double
dandc_sum( const double *a , 
	   const size_t lo ,
	   const size_t hi ) ;

/**
   @fn double kahan_summation( const double *a , const size_t LENGTH )
   @brief kahan summation algorithm
   @param a :: array of values
   @param LENGTH :: length of array a
   @return the sum of the values in @a
 */
double
kahan_summation( const double *a ,
		 const size_t LENGTH ) ;

/**
   @fn double knuth_average( const double *a , const size_t LENGTH )
   @brief (more) stable average attributed to knuth
   @param a :: array of values
   @param LENGTH :: length of array a
   @return the average of the values stored in @a
 */
double
knuth_average( const double *a , 
	       const size_t LENGTH ) ;

/**
   @fn double par_dandc_sum( const double *a , const size_t N )
   @brief recursive parralelised divide and conquer sum
   @param a :: array of values
   @param N :: length of array a
   @return the sum of the values in @a
   Note:: Error should go like Nlog(N)
 */
double
par_dandc_sum( const double *a , 
	       const size_t N ) ;

#endif
