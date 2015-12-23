/**
   @file GLU_sums.c
   @brief round-off resistant summations and averages
 */
#include <stdlib.h>

#ifdef HAVE_OMP_H
#include <omp.h>
#endif

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

// computes the variance as well as the average
void
knuth_average_err( double *ave , 
		   double *err , 
		   const double *data , 
		   const size_t N )
{
  size_t i ;
  *err = *ave = 0.0 ;
  for( i = 0 ; i < N ; i++ ) {
    const double delta = data[i] - *ave ;
    *ave += delta / ( i + 1.0 ) ;
    *err += delta * ( data[i] - *ave ) ;
  }
  *err /= ( N - 1 ) ;
  return ;
}

// parallel divide and conquer sum
#ifdef HAVE_OMP_H
double
par_danc_sum( const double *a , 
	      const size_t N )
{
  size_t i , Nth = 1 ;
  // open a parallel environment to get number of threads
  #pragma omp parallel
  {
    Nth = omp_get_num_threads() ;
  }
  const size_t split = ( N + Nth - 1 ) / Nth ;
  register double sum = 0 ;
  // anything less than the division goes to the last thread
  #pragma omp parallel for private(i) reduction(+:sum)
  for( i = 0 ; i < Nth ; i++ ) {
    sum = sum + (double)dandc_sum( a , i * split , 
				   (i+1)*split> N ? N : (i+1)*split ) ;
  }
  return sum ;
}
#endif
