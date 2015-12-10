/**
   @file LU.h
   @brief LU decomposition determinant
 */
#ifndef GLU_LU_H
#define GLU_LU_H

#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC ) && (NC>3)

// include the SSE version
#include "LU_SSE.h"

#else

/**
   @fn double complex LU_det( const int N , const GLU_complex U[ N*N ] )
   @brief computes the determinant of a matrix U
   @param N :: Side length of the square matrix
   @param U :: Matrix having its determinant taken
   @return the determinant
 */
double complex
LU_det( const int N , 
	const GLU_complex U[ N*N ] ) ;

/**
   @fn double complex LU_det_overwrite( const int N , GLU_complex U[ N*N ] )
   @brief computes the determinant of a matrix U
   @param N :: Side length of the square matrix
   @param U :: Matrix having its determinant taken, overwritten
   @return the determinant
   @warning overwrites space in U
 */
double complex
LU_det_overwrite( const int N , 
		  GLU_complex U[ N*N ] ) ;


#endif

#endif
