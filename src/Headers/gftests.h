/*
    Copyright 2013 Renwick James Hudspith

    This file (gftests.h) is part of GLU.

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
   @file gftests.h
   @brief function definitions for the gauge-fixing tests
 */

#ifndef GLU_GFTESTS_H
#define GLU_GFTESTS_H

/**
   @fn void const_time( const struct site *__restrict lat , double *const_lin , double *const_lin )
   @brief Constance in any direction suggests good landau gauge fixing
   @param lat :: gauge field
   @param const_lin :: Hermitian projection def of gauge fields t-constance
   @param const_log :: Exact log def of gauge fields t-constance

   @return returns the AntiHermitian_proj definition of the difference in the constance in time
 */
void
const_time( const struct site *__restrict lat ,
	    double *const_lin ,
	    double *const_log ) ;

/**
   @fn double gauge_functional( const struct site *__restrict lat )
   @brief computes the gauge fixing functional
   @param lat :: gauge fixed lattice links
   @return the functional value

   computes <\br>
   \f[
   1.0 - \frac{1}{ND*NC*LVOLUME} \sum_{x,\mu} \Re\left(tr\left[U_{\mu}(x)\right]\right)
   \f]
   for the hermitian projection definition of the links, and <\br>
   \f[
   \frac{1}{ND*NC*LVOLUME} \sum_{x,\mu} \Re\left(tr\left[A_\mu(x)^2\right]\right)
   \f]
   for the log-definition of the links.
 */
double
gauge_functional( const struct site *__restrict lat ) ;

/**
   @fn double gauge_test( const GLU_complex *__restrict *__restrict gauge )
   @brief Checks how close to unity our gauge transformation matrices are

   @param gauge :: Lattice-wide gauge transformation matrices

   The rationale here is that the gauge transformation matrices should be 1 when
   the gauge is completely fixed. Returns the lattice average of

   @return 
   \f[

   \frac{1}{NC \times VOLUME }\sum_x Tr( 1 - gauge(x) )

   \f]
 */
double
gauge_test( const GLU_complex *__restrict *__restrict gauge ) ;

/**
   @fn double gtrans_functional( const struct site *__restrict lat , const GLU_complex *__restrict *__restrict slice_gauge , const size_t t )
   @brief computes the SPATIAL functional for a given slice's gauge transformations
   @param lat :: lattice gauge fields
   @param slice_gauge :: slice-wide gauge transformation matrices
   @param t :: timeslice index
 */
double
gtrans_functional( const struct site *__restrict lat ,
		   const GLU_complex *__restrict *__restrict slice_gauge ,
		   const size_t t ) ;

/**
   @fn double theta_test_lin( const struct site *__restrict lat , GLU_real *max , const size_t MAX_DIR )
   @brief Checks the derivative of the AntiHermitian_proj definition of the lie fields
   
   @param lat :: gauge fields
   @param max :: maximum value of the derivative
   @param MAX_DIR :: maximum number of dimensions i.e. ND-1 would check
   configuration space coulomb gauge fixing criterion.

   @return
   \f[

   \frac{2}{NC\times VOLUME} \sum_x Tr | \partial_\mu A^{lin}_\mu (x) | ^2

   \f]
 */
double 
theta_test_lin( const struct site *__restrict lat , 
		GLU_real *max ,
		const size_t MAX_DIR ) ;

/**
   @fn double theta_test_log( const struct site *__restrict lat , GLU_real *max , const size_t MAX_DIR )
   @brief Checks the derivative of the exact log definition of the lie fields
   
   @param lat :: gauge fields
   @param max :: maximum value of the derivative
   @param MAX_DIR :: maximum number of dimensions i.e. ND-1 would check
   configuration space coulomb gauge fixing criterion.

   @return
   \f[

   \frac{2}{NC\times VOLUME} \sum_x Tr | \partial_\mu A^{log}_\mu (x) | ^2

   \f]
 */
double 
theta_test_log( const struct site *__restrict lat , 
		GLU_real *max ,
		const size_t MAX_DIR ) ;

#endif
