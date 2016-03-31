/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (lie_mats.h) is part of GLU.

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
   @file lie_mats.h
   @brief prototype functions for the computation of lie elements and structure functions for SU(NC) matrices
 */
#ifndef GLU_LIE_MATS_H
#define GLU_LIE_MATS_H

/**
   @fn void compute_fs_and_ds( void )
   @brief allocates the memory for, and computes the structure functions for fundamental SU(N) matrices
   
   Computes the structure functions using the formulae
   \f[
   f^{abc} = -2I TR\left[ [T^{a},T^{b}]T^{c}\right]   
   \f]
   \f[
   d^{abc} = 2 TR\left[ \{T^{a},T^{b}\} T^{c}\right]   
   \f]
   @warning allocates memory to static arrays, requires the initialisation of the generators init_generators(). Should be finished with a call to free_f_and_d(). 
 */
void
compute_fs_and_ds( void ) ;

/**
   @fn double complex dabc_ABC( const GLU_complex A[ NCNC ] , const GLU_complex B[ NCNC ] , const GLU_complex C[ NCNC ] )
   @brief contraction with the combination \f$ A^a B^b C^c (d^{abc}) \f$
   @param A :: Matrix A
   @param B :: Matrix B
   @param C :: Matrix C

   @warning matrices must be either hermitian or the Fourier transform of a hermitian matrix for this to make any sense at all.
   @return \f$ A^a B^b C^c (d^{abc}) \f$
 */
double complex
dabc_ABC( const GLU_complex A[ NCNC ] , 
	  const GLU_complex B[ NCNC ] , 
	  const GLU_complex C[ NCNC ] ) ;

/**
   @fn void free_generators( void )
   @brief frees the allocated generator matrices
 */
void
free_generators( void ) ;

/**
   @fn void free_f_and_d( void )
   @brief frees the memory of the computed f and d structure functions
 */
void
free_f_and_d( void ) ;

/**
   @fn double complex ifabc_ABC( const GLU_complex A[ NCNC ] , const GLU_complex B[ NCNC ] , const GLU_complex C[ NCNC ] )
   @brief contraction with the combination \f$ A^a B^b C^c (if^{abc}) \f$
   @param A :: Matrix A
   @param B :: Matrix B
   @param C :: Matrix C

   @warning matrices must be either hermitian or the Fourier transform of a hermitian matrix for this to make any sense at all.
   @return \f$ A^a B^b C^c (if^{abc}) \f$
 */
double complex
ifabc_ABC( const GLU_complex A[ NCNC ] , 
	  const GLU_complex B[ NCNC ] , 
	  const GLU_complex C[ NCNC ] ) ;

/**
   @fn double complex ifabc_dabc_ABC( const GLU_complex A[ NCNC ] , const GLU_complex B[ NCNC ] , const GLU_complex C[ NCNC ] )
   @brief contraction with the combination \f$ A^a B^b C^c (I.f^{abc}+d^{abc}) 4.Tr(ABC) \f$
   @param A :: Matrix A
   @param B :: Matrix B
   @param C :: Matrix C

   @warning matrices must be either hermitian or the Fourier transform of a hermitian matrix for this to make any sense at all.
   @return \f$ A^a B^b C^c (I.f^{abc}+d^{abc}) \f$
 */
double complex
ifabc_dabc_ABC( const GLU_complex A[ NCNC ] , 
		const GLU_complex B[ NCNC ] , 
		const GLU_complex C[ NCNC ] ) ;

/**
   @fn void init_generators( void )
   @brief initialises the generator matrices for fundamental SU(N) matrices
   @warning allocates memory, must be used in conjunction with free_generators()
 */
void
init_generators( void ) ;

#endif
