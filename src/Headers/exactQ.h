/*
Copyright 2013-2025 Renwick James Hudspith

    This file (exactQ.h) is part of GLU.

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
   @file exactQ.h
   @brief Function definitions for calculating the lie-matrices from links.
 */
#ifndef GLU_EXACTQ_H
#define GLU_EXACTQ_H

/**
   @fn void AntiHermitian_proj( GLU_complex Q[ NCNC ] , const GLU_complex U[ NCNC ] )
   @brief Projection of an 3x3 matrix to its traceful-antihermitian part...
   @param Q [out] :: The lie matrix calculated
   @param U [in] :: The link matrix
   @warning overwrites Q
   
   Performs 
   \f[ 
   Q(x) = \frac{1}{2}\left[\left( U(x)-U(x)^{\dagger}\right)\right]
   \f]
 **/
void 
AntiHermitian_proj( GLU_complex Q[ NCNC ] , 
		    const GLU_complex U[ NCNC ] ) ;

/**
   @fn void AntiHermitian_proj_short( GLU_complex Q[ HERMSIZE ] , const GLU_complex U[ NCNC ] )
  @brief Projection of an 3x3 matrix to its traceful-antihermitian short part...
  @param Q [out] :: The lie matrix calculated
  @param U [in] :: The link matrix
  @warning overwrites Q

   Performs 
   \f[ 
   Q(x) = \frac{1}{2}\left[\left( U(x)-U(x)^{\dagger}\right)\right]
   \f]

  Returns the shortened representation of the hermitian
  matrix i.e. the upper diagonal as the bottom is related through conjugacy
 **/
void 
AntiHermitian_proj_short( GLU_complex Q[ HERMSIZE ] , 
	      const GLU_complex U[ NCNC ] ) ;

/**
   @fn void exact_log_slow( GLU_complex Q[ NCNC ] , const GLU_complex U[ NCNC ] )
   @brief Slow exact logarithm
   @param Q [out] :: The exact lie-field
   @param U [in] :: The field we are taking the exact logarithm of
   @warning overwrites a
   uses exact exponentiation methods in effs.h
 **/
void 
exact_log_slow( GLU_complex Q[ NCNC ] , 
		const GLU_complex U[ NCNC ] ) ;

/**
    @fn void exact_log_slow_short( GLU_complex Q[ HERMSIZE ] , const GLU_complex U[ NCNC ] )
   @brief Slow exact logarithm shortened output
   @param Q :: The exact lie-field
   @param U :: The field we are taking the exact logarithm of
   @warning overwrites a
   Returns the upper diagonal of the matrix. uses exact exponentiation methods in effs.h
 **/
void 
exact_log_slow_short( GLU_complex Q[ HERMSIZE ] , 
		      const GLU_complex U[ NCNC ] ) ;

/**
   @fn void get_iQ( GLU_complex Q[ NCNC ] , const GLU_complex U[ NCNC ] )
   @brief Uses "v^{\dagger}U v" to diagonalise U. Takes log of this and multiplies through again by "v, v^{\dagger}" to get the log Q
   @param U :: Matrix we are taking the log of
   @param Q :: The log of U
 */
void 
get_iQ( GLU_complex Q[ NCNC ] ,
	const GLU_complex U[ NCNC ] ) ;

/**
   @fn void Hermitian_proj( GLU_complex Q[ NCNC ] , const GLU_complex U[ NCNC ] )
   @brief Projection of an 3x3 matrix to its traceless hermitian part...
   @param Q [out] :: The lie matrix calculated
   @param U [in] :: The link matrix
   
   @warning overwrites Q
   
   Performs 
   \f[ 
   Q(x) = \frac{1}{2i}\left[\left( U(x)-U(x)^{\dagger}\right) - \frac{1}{NC}\left( U(x)-U(x)^{\dagger}\right)I_{NC \times NC}\right]
   \f]
**/
void 
Hermitian_proj( GLU_complex Q[ NCNC ] , 
		const GLU_complex U[ NCNC ] ) ;

/**
   @fn void Hermitian_proj_short( GLU_complex Q[ HERMSIZE ] , const GLU_complex U[ NCNC ] )
  @brief Projection of an 3x3 matrix to its traceless hermitian part...
  @param Q [out] :: The lie matrix calculated
  @param U [in] :: The link matrix

  @warning overwrites Q

  Performs 
  \f[ 
  Q(x) = \frac{1}{2i}\left[\left( U(x)-U(x)^{\dagger}\right) - \frac{1}{NC}\left( U(x)-U(x)^{\dagger}\right)I_{NC \times NC}\right]
  \f]

  Returns the shortened representation of the hermitian
  matrix i.e. the upper diagonal as the bottom is related through conjugacy
 **/
void 
Hermitian_proj_short( GLU_complex Q[ HERMSIZE ] , 
		      const GLU_complex U[ NCNC ] ) ;

/**
   @fn void trf_AntiHermitian_proj( GLU_complex Q[ NCNC ] , const GLU_complex U[ NCNC ] )
  @brief Projection of an 3x3 matrix to its traceful-antihermitian part...
  @param Q :: The lie matrix calculated
  @param U :: The link matrix
  @warning overwrites Q

  Performs 
  \f[ 
  Q(x) = \frac{1}{2}\left[\left( U(x)-U(x)^{\dagger}\right) - \frac{1}{NC}\left( U(x)-U(x)^{\dagger}\right)I_{NC \times NC}\right]
  \f]
 **/
void 
trf_AntiHermitian_proj( GLU_complex Q[ NCNC ] , 
			const GLU_complex U[ NCNC ] ) ;

#endif
