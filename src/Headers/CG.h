/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (CG.h) is part of GLU.

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
   @file CG.h
   @brief common routines between Landau and Coulomb for the CG routines
 */
#ifndef GLU_CG_H
#define GLU_CG_H

/**
   @fn void allocate_traces( const size_t LENGTH )
   @brief allocate the traces array that is used in the average
   @param LENGTH :: the length of the traces array you want to allocate
 */
void allocate_traces( const size_t LENGTH ) ;

/**
   @fn double approx_minimum( const size_t nmeas , const double alphas[ nmeas ] , const double functional[ nmeas ] )
   @brief finds the approximate minimum of alphas using GLU-bic splines
   @param nmeas :: number of measurements made
   @param alphas :: alphas tested
   @param functional :: the gauge functional at each alpha

   @return the alpha that approximately minimises the functional, or 0 
 */
double
approx_minimum( const size_t nmeas , 
		const double alphas[ nmeas ] ,
		const double functional[ nmeas ] ) ;

/**
   @fn double coul_gtrans_fields( struct sp_site_herm *__restrict rotato , const struct site *__restrict lat , const GLU_complex *__restrict *__restrict slice_gauge , const size_t t )
   @brief computes the gauge transformed Lie fields
   @param rotato :: (Coulomb) gauge rotated Lie fields
   @param lat :: link matrices
   @param slice_gauge :: gauge transformation matrices on this time-slice
   @param t :: time-slice index
   @param acc :: gauge fixing accuracy
 */
double
coul_gtrans_fields( struct sp_site_herm *__restrict rotato ,
		    const struct site *__restrict lat ,
		    const GLU_complex *__restrict *__restrict slice_gauge ,
		    const size_t t ) ;

/**
   @fn double evaluate_alpha( const GLU_complex *__restrict *__restrict gauge ,	const struct site *__restrict lat , const size_t DIR ,	const size_t LENGTH , const size_t t ) 
   @brief evaluate the functional for a test expansion parameter alpha
   @param gauge :: gauge transformation matrices
   @param lat :: lattice links
   @param DIR :: (#ND-1) Coulomb or (#ND) for Landau
   @param LENGTH :: the length of the gauge transformation array
   @param t :: the timeslice index (set to 0 for Landau)
   @return the average functional
 */
double
evaluate_alpha( const GLU_complex *__restrict *__restrict gauge , 
		const struct site *__restrict lat ,
		const size_t DIR ,
		const size_t LENGTH ,
		const size_t t ) ;

/**
   @fn void FOURIER_ACCELERATE( GLU_complex *__restrict *__restrict in , GLU_complex *__restrict *__restrict out , const void *__restrict forward , const void *__restrict backward , const double *__restrict psq , const size_t LENGTH )
   @brief perform a Fourier Acceleration correction to the derivative in
   @param in :: derivative of the lie fields
   @param out :: FFT storage
   @param forward :: FFTW forward plan
   @param backward :: FFTW backward plan
   @param psq :: is \f$ i \frac{ p^{2}_{Max} }{LENGTH p^2} \f$
   @param LENGTH :: is the LENGTH of all the arrays #LCU for Coulomb #LVOLUME for Landau
 */
void
FOURIER_ACCELERATE( GLU_complex *__restrict *__restrict in ,
		    GLU_complex *__restrict *__restrict out ,
		    const void *__restrict forward ,
		    const void *__restrict backward ,
		    const GLU_real *__restrict psq ,
		    const size_t LENGTH ) ;

/**
   @fn free_traces( )
   @brief deallocate the static "traces" array in CG.c
 */
void free_traces( void ) ;

/**
   @fn double gauge_functional_fast( const struct site *__restrict lat )
   @brief performs a fast evaluation of the functional
   @param lat :: the lattice links
   Assumes lat has been gauge transformed, this is used in the alpha=0 case
 */
double
gauge_functional_fast( const struct site *__restrict lat ) ;

/**
   @fn double PRfmax( const double a , const double b )
   @brief maximum of a and b
 */
double
PRfmax( const double a , 
	const double b ) ;

/**
   @fn void set_gauge_matrix( GLU_complex *__restrict gx , const GLU_complex *__restrict *__restrict in , const double alpha , const size_t i ) 
   @brief unpacks and exponentiates the derivative into the array gx at site i
   @param gx :: gauge transformation matrix
   @param in :: derivative at site i
   @param i :: site index
 */
void
set_gauge_matrix( GLU_complex *__restrict gx ,
		  const GLU_complex *__restrict *__restrict in ,
		  const double alpha ,
		  const size_t i ) ;

/**
   @fn double sum_deriv( const GLU_complex *__restrict *__restrict in , const size_t LENGTH )
   @brief computes the sum of the trace of the product of the derivative matrices
   @param in :: derivative matrix
   @param LENGTH :: length of the array in
   
   @return \f$ \sum_x Tr\left[ (\Delta_{\mu} A_{\mu})^2 \right] \f$
 */
double
sum_deriv( const GLU_complex *__restrict *__restrict in , 
	   const size_t LENGTH ) ;

/**
   @fn double sum_PR_numerator( const GLU_complex *__restrict *__restrict in , const GLU_complex *__restrict *__restrict in_old , const size_t LENGTH )
   @brief the sum in * ( in - in_old ), is the Polak-Ribiere numerator
   @param in :: the derivative at this step
   @param in_old :: the previous derivative
   @param LENGTH :: length of the arrays
   @warning in and in_old should be #TRUE_HERM length arrays in FFTW order
 */
double
sum_PR_numerator( const GLU_complex *__restrict *__restrict in , 
		  const GLU_complex *__restrict *__restrict in_old ,
		  const size_t LENGTH ) ;

#endif
