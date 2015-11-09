/**
   @file GLU_sums.h
   @brief prototype declarations for robust summations
 */
#ifndef GLU_SUMS_H
#define GLU_SUMS_H

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
