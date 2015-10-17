/**
   @file GLU_sums.c
   @brief round-off resistant summations and averages
 */
#include <stdlib.h>

// numerically more friendly average?
double
knuth_average( const double *a , 
	       const size_t LENGTH )
{
  size_t i ;
  register double ave = 0.0 ;
  for( i = 0 ; i < LENGTH ; i++ ) {
    ave += ( a[i] - ave ) / ( 1.0 + i ) ;
  }
  return ave ;
}

// kahan summation
double
kahan_summation( const double *a ,
		 const size_t LENGTH )
{
  size_t i ;
  register double sum = 0.0 , c = 0.0 , y , t ;
  for( i = 0 ; i < LENGTH ; i++ ) {
    y = a[i] - c ;
    t = sum + y ;
    c = ( t - sum ) - y ;
    sum = t ;
  }
  return sum ;
}

// divide and conquer summation
double
dandc_sum( const double *a , 
	   const size_t lo ,
	   const size_t hi )
{
  return ( hi - lo ) < 2 ? a[lo] : \
    dandc_sum( a , lo , (hi+lo)/2 ) +\
    dandc_sum( a , (hi+lo)/2 , hi ) ;
}
