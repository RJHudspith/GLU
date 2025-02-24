/*
Copyright 2013-2025 Renwick James Hudspith

    This file (effs.h) is part of GLU.

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
   @file effs.h
   @brief Function definitions for the exact exponentiation of lie fields
 */
#ifndef GLU_EFFS_H
#define GLU_EFFS_H

/**
   @fn void calculate_effs_VDM_herm( double complex *__restrict f , const double *__restrict z )
   @brief computes the f-constants from the general Vandermonde system
   @param f [out] :: f-constants
   @param z [in] :: Eigenvalues 

   Takes the eigenvalues as input and returns the f's. Adapted from Golub and Van Loan
 */
void 
calculate_effs_VDM_herm( double complex *__restrict f , 
			 const double *__restrict z ) ;

/**
   @fn void calculate_effs_VDM_suNC( double complex *__restrict f , const double complex *__restrict z )
   @brief computes the f-constants from the general Vandermonde system
   @param f [out] :: f-constants
   @param z [in] :: Eigenvalues 

   Takes the eigenvalues as input and returns the f's. Adapted from Golub and Van Loan
 */
void 
calculate_effs_VDM_suNC( double complex *__restrict f , 
			 const double complex *__restrict z ) ;


/**
   @fn void f_hermitian_log_suNC( double complex f[NC] , const double complex z[NC] )
   @param f :: The fs
   @param z :: The eigenvalues of the U matrix

   SU3 only. Function for taking the log of an SU(#NC) matrix by calculating its
   f's which are the same as the lie-matrix's.

   @warning I provide a cut off at the ulp of GLU_real precision to limit numerical instability

   @warning SU(NC) ONLY!
 */
#if NC == 2
/**
   @fn void f_hermitian_log_suNC( double complex f[NC] , const double z )
   @param f :: The fs
   @param z :: The eigenvalue of the U matrix, su2 so +/-
   @warning I provide a cut off at the ulp of GLU_real precision to limit numerical instability
void 
*/
void
f_hermitian_log_suNC( double complex f[NC] , 
		      const double z ) ;
#else
/**
   @fn void f_hermitian_log_suNC( double complex f[NC] , const double complex z[NC] )
   @param f :: The fs
   @param z :: The eigenvalues of the U matrix
   @warning I provide a cut off at the ulp of GLU_real precision to limit numerical instability
 */
void 
f_hermitian_log_suNC( double complex f[NC] , 
		      const double complex z[NC] ) ;
#endif

#endif
