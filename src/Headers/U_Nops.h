/*
    Copyright 2013 Renwick James Hudspith

    This file (U_Nops.h) is part of GLU.

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
   @file U_Nops.h
   @brief matrix operations live here
 */
#ifndef GLU_U_NOPS_H
#define GLU_U_NOPS_H

/**
   @fn void add_constant( GLU_complex a[ NCNC ] , const GLU_complex c ) ;
   @brief atomically adds to a the constant on the diagonal "c"
   @param a :: a square matrix
   @param c :: generic constant
   @warning atomic addition into @a
 */
void
add_constant( GLU_complex a[ NCNC ] , 
	      const GLU_complex c ) ;

/**
   @fn void a_plus_b( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] )
   @brief atomic matrix addition into a, e.g a[i] += b[i] all i
   @param a :: matrix having "@a b" added
   @param b :: matrix being added to @a a
   @warning atomically overwrites @a a
 */
void
a_plus_b( GLU_complex a[ NCNC ] , 
	  const GLU_complex b[ NCNC ] ) ;
/**
   @fn void a_plus_CSxb( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex S )
   @brief computes a += S * b , where S is some GLU_complex -valued scalar
   @param a :: matrix having elements atomically added
   @param b :: arbitrary NC x NC matrix
   @param S :: GLU_complex valued scalar
   @warning S has to be a GLU_complex!
 */
void
a_plus_CSxb( GLU_complex a[ NCNC ] , 
	  const GLU_complex b[ NCNC ] ,
	  const GLU_complex S ) ;

/**
   @fn void a_plus_Sxb( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_real S )
   @brief computes a += S * b , where S is some GLU_real-valued scalar
   @param a :: matrix having elements atomically added
   @param b :: arbitrary NC x NC matrix
   @param S :: GLU_real valued scalar
   @warning @a S has to be a GLU_real! @a a is atomically added!
 **/
void
a_plus_Sxb( GLU_complex a[ NCNC ] , 
	  const GLU_complex b[ NCNC ] ,
	  const GLU_real S ) ;

/**
   @fn void a_plus_Sxbminc_short( GLU_complex a[ HERMSIZE ] , const GLU_real S , const GLU_complex B[ HERMSIZE ] , const GLU_complex C[ HERMSIZE ] )
   @brief computes the shortened a = S * ( b - c ) where a,b and c are matrices
   @param a :: result of the operation
   @param S :: real constant
   @param B :: some #HERMSIZE matrix
   @param C :: some other #HERMSIZE matrix
 */
void
a_plus_Sxbminc_short( GLU_complex a[ HERMSIZE ] ,
		      const GLU_real S ,
		      const GLU_complex B[ HERMSIZE ] ,
		      const GLU_complex C[ HERMSIZE ] ) ;

/**
   @fn void b_min_c( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] ,const GLU_complex c[ NCNC ] )
   @param a :: output matrix
   @param b :: arbitrary NC x NC GLU_complex matrix
   @param c :: arbitrary NC x NC GLU_complex matrix
   @brief computes matrix negation \f$ a = b - c \f$ , element by element
   @warning overwrites a
 */
void 
b_min_c( GLU_complex a[ NCNC ] , 
	 const GLU_complex b[ NCNC ] ,
	 const GLU_complex c[ NCNC ] ) ;

/**
   @fn void cofactor_transpose( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] )
   @brief computes the cofactor transpose of matrix b and writes into a
   @param a :: transpose of the cofactor matrix
   @param b :: matrix having its cofactor transpose taken
   @warning a is overwritten
   the cofactor transpose transposes the cofactor matrix. Each element of the cofactor matrix is built from the signed minors of not in the column and row of that element

   @return the determinant
 */
GLU_complex
cofactor_transpose( GLU_complex a[ NCNC ] ,
		    const GLU_complex b[ NCNC ] ) ;

/**
   @fn void equiv( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] )
   @brief equates matrix a to b
   @param a :: is overwritten with b
   @param b :: what we are equating a with
   @warning they should be the same size
 */
void 
equiv( GLU_complex a[ NCNC ] , 
       const GLU_complex b[ NCNC ] ) ;

/**
   @fn void dagger( GLU_complex b[ NCNC ] , const GLU_complex a[ NCNC ] )
   @brief computes the conjugate transpose \f$ a = b^{\dagger} \f$
   @param a :: matrix outputted
   @param b :; matrix we are conjugate-transposing
   @warning overwrites a
 */
void 
dagger( GLU_complex b[ NCNC ] , 
	const GLU_complex a[ NCNC ] ) ;

/**
   @fn void diag( GLU_complex M[ NCNC ] , const GLU_complex c ) 
   @brief turns "M" into a diagonal matrix with diagonals of constant "c"
   @param M :: diagonal matrix of constant c
   @param c :: constant
   @warning overwrites M
   off diagonal elements are zeroed
 */
void
diag( GLU_complex M[ NCNC ] ,
      const GLU_complex c ) ;

/**
   @fn void diag_vect( GLU_complex M[ NCNC ] , const GLU_complex c[ NC ] ) 
   @brief turns "M" into a diagonal matrix with diagonals of constant vector "c"
   @param M :: diagonal matrix of constant vector c[ NC ]
   @param c :: constant vector [ NC ] 
   @warning overwrites M
   off diagonal elements are zeroed <br>
   for a 3x3 matrix M[0] = c[0], M[4] = c[1] , M[8] = c[2] rest 0
 */
void
diag_vect( GLU_complex M[ NCNC ] ,
	   const GLU_complex c[ NC ] ) ;

/**
   @fn void identity( GLU_complex ident[ NCNC ] )
   @brief overwrites the matrix ident to be the identity
   @param ident :: matrix being set to the identity
   @warning overwrites matrix ident
 */
void 
identity( GLU_complex ident[ NCNC ] ) ;

/**
   @fn int is_unitary( const GLU_complex U[ NCNC ] )
   @brief check if the matrix is unitary
   @param U :: matrix being checked
   @return #GLU_SUCCESS or #GLU_FAILURE
   @warning checks to a precision of \f$ 10^{-6} \f$ for Single prec and \f$ 10^{-14} \f$ for double.
 */
GLU_bool
is_unitary( const GLU_complex U[ NCNC ] ) ;

/**
   @fn GLU_complex det( const GLU_complex U[ NCNC ] )
   @brief computes the determinant of a matrix U
   @param U :: Matrix having its determinant taken
   @return the determinant
   does not demand any symmetry to work, good test of SU(N)-ness
 */
GLU_complex 
det( const GLU_complex U[ NCNC ] ) ;

/**
   @fn void mat_mult_vec( GLU_complex vect[ NC ] , const GLU_complex S[ NCNC ] , const GLU_complex v[ NC ] )
   @brief matrix multiplied by a vector, e.g. vect = S·v
   @param vect :: output vector
   @param S :: arbitrary NC x NC matrix
   @param v :: vector being multiplied by @a S
   @warning should be all the same size
 */
void 
mat_mult_vec( GLU_complex vect[ NC ] , 
	      const GLU_complex S[ NCNC ] , 
	      const GLU_complex v[ NC ] ) ;

/**
   @fn void M_times_c( GLU_complex M[ NCNC ] , const GLU_complex c )
   @brief multiplies matrix "M" by constant "c"
   @param M :: matrix being multiplied by c
   @param c :: arbitrary GLU_complex constant
   @warning overwrites the matrix "M"
 */
void 
M_times_c( GLU_complex M[ NCNC ] , 
	   const GLU_complex c ) ;

/**
   @fn void matrix_power( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const int n )
   @brief fastpower-like matrix multiply routine
   @param a :: the result a = b^n
   @param b :: the matrix being raised to a power n
   @param n :: the integer power
   @warnng only works for non-negative integer powers
 */
void
matrix_power( GLU_complex a[ NCNC ] , 
	      const GLU_complex b[ NCNC ] , 
	      const int n ) ;

/**
   @fn void outerproduct( GLU_complex Q[ NCNC ] , const GLU_complex a[ NC ] , const GLU_complex b[ NC ] )
   @brief computes the outerproduct of two vectors \f$ Q = a \otimes b \f$
   @param Q :: the outerproduct of a and b
   @param a :: arbitrary vector
   @param b :: arbitrary vector 
   @warning a and be should be of length NC
 */
void 
outerproduct( GLU_complex Q[ NCNC ] , 
	      const GLU_complex a[ NC ] , 
	      const GLU_complex b[ NC ] ) ;

/**
   @fn void pack_hermitian( GLU_complex a[ HERMSIZE ] , const GLU_complex b[ NCNC ] ) 
   @brief takes a Hermitian matrix b and packs it into the HERMSIZE a 
   @param a :: packed matrix #HERMSIZE
   @param b :: hermitian matrix
 */
void 
pack_hermitian( GLU_complex a[ HERMSIZE ] , 
		const GLU_complex b[ NCNC ] ) ;

/**
   @fn void printcomplex( const GLU_complex a )
   @brief prints to stdout the complex variable "a"
   @param a :: prints this out
 */
void 
printcomplex( const GLU_complex a ) ;

/**
   @fn void rebuild( GLU_complex a[ NCNC ] , const GLU_real *__restrict b )
   @brief rebuilds the shortened matrix form of a into its full form
   @param a :: SU(N) matrix
   @param b :: parameters the matrix has been shortened to
 */
void 
rebuild( GLU_complex a[ NCNC ] , 
	 const GLU_real *__restrict b ) ;

/**
   @fn void rebuild_antihermitian( GLU_complex a[ NCNC ] , const GLU_complex b[ HERMSIZE ] )
   @brief takes the upper-diagonal entries "b" and conjugates them into the lower triangular of "a" and puts "b" into the upper triangular
   @param a :: full antihermitian matrix
   @param b :: our #HERMSIZE representation
   @warning tracelessness enforced
 */
void 
rebuild_antihermitian( GLU_complex a[ NCNC ] , 
		       const GLU_complex b[ HERMSIZE ] ) ;

/**
   @fn void rebuild_hermitian( GLU_complex a[ NCNC ] , const GLU_complex b[ HERMSIZE ] )
   @brief takes the upper-diagonal entries "b" and conjugates them into the lower triangular of "a" and puts "b" into the upper triangular
   @param a :: full hermitian matrix
   @param b :: our #HERMSIZE representation
   @warning tracelessness enforced
 */
void 
rebuild_hermitian( GLU_complex a[ NCNC ] , 
		   const GLU_complex b[ HERMSIZE ] ) ;

/**
   @fn void shorten( GLU_real *__restrict a , const GLU_complex b[ NCNC ] )
   @brief shortens the special-unitary matrix b to its smallest representation
   @param a :: parameters the matrix has been shortened to
   @param b :: SU(N) matrix
   @warning does not work if b is exactly unity
   for SU(3) a is 8 real parameters and for SU(2) it is 3 real parameters
 */
void
shorten( GLU_real *__restrict a ,
	 const GLU_complex b[ NCNC ] ) ;

/**
   @fn void speed_det( GLU_complex *__restrict dt ,  const GLU_complex U[ NCNC ] )
   @brief computes the determinant of matrix "U"
   @param U :: Matrix having its determinant taken
   @param dt :: the result passed by reference
 */
void 
speed_det( GLU_complex *__restrict dt ,
	   const GLU_complex U[ NCNC ] ) ;

/**
   @fn void speed_trace( GLU_complex *res , const GLU_complex U[ NCNC ] )
   @brief computes the trace of "U" and passes by reference the result
   @param U :: the matrix we are tracing
   @param res :: the result
 */
void 
speed_trace( GLU_complex *res ,
	     const GLU_complex U[ NCNC ]  ) ;

/**
   @fn void speed_trace_Re( double *__restrict res , const GLU_complex U[ NCNC ] )
   @brief computes the trace of "U" and passes by reference the real part of the result 
   @param U :: the matrix we are tracing
   @param res :: the result
 */
void 
speed_trace_Re( double *__restrict res ,
		const GLU_complex U[ NCNC ] ) ;

/**
   @fn GLU_complex trace( const GLU_complex U[ NCNC ] )
   @brief computes the trace of matrix "U"
   @param U :: matrix we are taking the trace of
   @return the trace
 */
GLU_complex 
trace( const GLU_complex U[ NCNC ] ) ;

/**
   @fn void trace_ab( GLU_complex *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC] )
   @brief computes the trace of the product of two matrices tr[ ab ]
   @param tr :: result passed by reference
   @param a :: arbitrary NC x NC matrix
   @param b :: arbitrary NC x NC matrix

   does not care about the type of matrix

   @warning a and b should be NC-square matrices
 */
void
trace_ab( GLU_complex *__restrict tr , 
	  const GLU_complex a[ NCNC ] , 
	  const GLU_complex b[ NCNC] ) ;

/**
   @fn void trace_abc( GLU_complex *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC] , const GLU_complex c[ NCNC ] )
   @brief computes the trace of the product of three matrices \f$ tr[ abc^{\dagger} ] \f$
   @param tr :: result passed by reference
   @param a :: arbitrary NC x NC matrix
   @param b :: arbitrary NC x NC matrix
   @param c :: arbitrary NC x NC matrix

   does not care about the type of matrix. Used for computing the trace of a gauge transformation.

   @warning a,b and c should be NC-square matrices
 */
void
trace_abc( GLU_complex *__restrict tr , 
	   const GLU_complex a[ NCNC ] , 
	   const GLU_complex b[ NCNC ] , 
	   const GLU_complex c[ NCNC ] ) ;

/**
   @fn void trace_abc_dag( GLU_complex *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC] , const GLU_complex c[ NCNC ] )
   @brief computes the trace of the product of three matrices tr[ abc ] where c is daggered
   @param tr :: result passed by reference
   @param a :: arbitrary NC x NC matrix
   @param b :: arbitrary NC x NC matrix
   @param c :: arbitrary NC x NC matrix

   does not care about the type of matrix

   @warning a,b and c should be NC-square matrices
 */
void
trace_abc_dag( GLU_complex *__restrict tr , 
	       const GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC ] , 
	       const GLU_complex c[ NCNC ] ) ;

/**
   @fn void trace_abc_dag_Re( GLU_real *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC] , const GLU_complex c[ NCNC ] )
   @brief computes the trace of the product of three matrices tr[ abc ] where c is daggered
   @param tr :: real part of the result passed by reference
   @param a :: arbitrary NC x NC matrix
   @param b :: arbitrary NC x NC matrix
   @param c :: arbitrary NC x NC matrix

   does not care about the type of matrix

   @warning a,b and c should be NC-square matrices
 */
void
trace_abc_dag_Re( GLU_real *__restrict tr , 
		  const GLU_complex a[ NCNC ] , 
		  const GLU_complex b[ NCNC ] , 
		  const GLU_complex c[ NCNC ] ) ;

/**
   @fn void trace_ab_dag( GLU_complex *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC] )
   @brief computes the trace of the product of two matrices \f$ tr\left[ ab^{\dagger} \right] \f$ where b is daggered
   @param tr :: result passed by reference
   @param a :: arbitrary NC x NC matrix
   @param b :: arbitrary NC x NC matrix

   does not care about the type of matrix

   @warning a and b should be NC-square matrices
 */
void
trace_ab_dag( GLU_complex *__restrict tr , 
	      const GLU_complex a[ NCNC ] , 
	      const GLU_complex b[ NCNC] ) ;

/**
   @fn void trace_ab_dag_Re( GLU_complex *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC] )
   @brief computes the trace of the product of two matrices \f$ tr\left[ ab^{\dagger} \right] \f$ where b is daggered
   @param tr :: real part of the result passed by reference
   @param a :: arbitrary NC x NC matrix
   @param b :: arbitrary NC x NC matrix

   does not care about the type of matrix

   @warning a and b should be NC-square matrices
 */
void
trace_ab_dag_Re( GLU_real *__restrict tr , 
		 const GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC] ) ;

/**
   @fn void trace_ab_herm( GLU_real *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC] )
   @brief trace of the product of two hermitian matrices, is always real
   @param tr :: the trace, passed by reference
   @param a :: matrix "a"
   @param b :: matrix "b"
 */
void
trace_ab_herm( GLU_real *__restrict tr , 
	       const GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC] ) ;

/**
   @fn void trace_ab_herm_short( GLU_real *__restrict tr , const GLU_complex a[ HERMSIZE ] , const GLU_complex b[ HERMSIZE ] )
   @brief computes the trace of the product of two hermitian short matrices
   @warning @a and @b have to Hermitian
 */
void
trace_ab_herm_short( GLU_real *__restrict tr , 
		     const GLU_complex a[ HERMSIZE ] , 
		     const GLU_complex b[ HERMSIZE ] ) ;

/**
   @fn void trace_prod_herm( GLU_real *__restrict tr , const GLU_complex a[ NCNC ] ) 
   @brief computes the trace of the product of two traceless hermitian matrices
   @param tr :: result passed by reference
   @param a :: hermitian matrix
   The answer has to be purely real as it is the sum of the absolute value of every element

   @warning @a has to be hermitian, hermiticity not checked
 */
void 
trace_prod_herm( GLU_real *__restrict tr ,
		 const GLU_complex a[ NCNC ] ) ;

/**
   @fn void transpose( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] )
   @brief computes the transpose \f$ a = b^{T} \f$
   @param a :: matrix outputted
   @param b :; matrix we are transposing
   @warning overwrites a
 */
void 
transpose( GLU_complex a[ NCNC ] , 
	   const GLU_complex b[ NCNC ] ) ;

/**
   @fn void write_matrix( const GLU_complex U[ NCNC ] )
   @brief writes to stdout the matrix "U"
   @param U :: matrix being written to the screen
   @warning prints out "%f" precision
 */
void 
write_matrix( const GLU_complex U[ NCNC ] ) ;

/**
   @fn void write_matrix_cform( const GLU_complex U[ NCNC ] )
   @brief writes to stdout the matrix "U" in the c GLU_complex array form
   @param U :: matrix being written to the screen
   @warning prints out "%f" precision
 */
void 
write_matrix_cform( const GLU_complex U[ NCNC ] ) ;

/**
   @fn void write_matrix_mathematica( const GLU_complex U[ NCNC ] )
   @brief writes out the matrix "U" to stdout in the format that can be copy-pasted into mathematica for operation-checking
   @param U :: matrix being written to the screen
   @warning prints out "%1.15f" precision
 */
void 
write_matrix_mathematica( const GLU_complex U[ NCNC ] ) ;

/**
   @fn void zero_mat( GLU_complex a[ NCNC ] )
   @brief sets the matrix "a" to zero
   @param a :: a square matrix
 */
void 
zero_mat( GLU_complex a[ NCNC ] ) ;

#endif
