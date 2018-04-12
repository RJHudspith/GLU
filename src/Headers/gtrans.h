/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (gtrans.h) is part of GLU.

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
   @file gtrans.h
   @brief function definitions for computing the gauge transformation and fourier accelerated landau gauge fixing routines.
 */

#ifndef GLU_GTRANS_H
#define GLU_GTRANS_H

/**
   @fn void gtransform_local( const GLU_complex *__restrict a , GLU_complex *__restrict b , const GLU_complex *__restrict c )
   @brief computes a.b.c^{\dagger}
 */
void
gtransform_local( const GLU_complex *__restrict a ,
		  GLU_complex *__restrict b ,
		  const GLU_complex *__restrict c ) ;

/**
   @fn void gtransform( struct site *__restrict lat , const GLU_complex *__restrict *__restrict gauge )
   @brief the lattice wide gauge transform. 

   @param lat :: lattice gauge fields
   @param gauge :: lattice gauge transformation matrices
   Taking the gauge matrices
   at the points and using them to rotate the link matrices.
   Implements the formula <br>

   @warning overwrites @a lat

   \f[
   U_{\mu}'(x) = g(x)U_{\mu}(x)g(x+\mu)^{\dagger}
   \f]
   \f[
   g(x),U_\mu(x)\in SU(NC)
   \f]
 **/
void 
gtransform( struct site *__restrict lat ,
	    const GLU_complex *__restrict *__restrict gauge ) ;

void 
gtransform2( struct site *__restrict lat ,
	     const GLU_complex *__restrict *__restrict gauge ) ;

/**
   @fn void gtransform_slice( const GLU_complex *__restrict *__restrict gauge , struct site *__restrict lat , const GLU_complex *__restrict *__restrict gauge_up , const size_t t )
   @brief The gauge transform routine used in the coulomb gauge fixing routine
   @param gauge :: gauge transformation matrices for this timeslice
   @param lat :: lattice gauge field
   @param gauge_up :: gauge transformation matrices for the timeslice above
   @param t :: timeslice index

   I perform the coulomb gauge fixing slice-by-slice, so I only need this type 
   of transform 

   @warning overwrites @a lat

   \f[
   U_{\mu}'(x) = g(x)U_{\mu}(x)g(x+\mu)^{\dagger}
   \f]
   \f[
   g(x),U_\mu(x)\in SU(NC)
   \f]
 **/
void
gtransform_slice( const GLU_complex *__restrict *__restrict gauge , 
		  struct site *__restrict lat , 
		  const GLU_complex *__restrict *__restrict gauge_up ,
		  const size_t t ) ;

/**
   @fn void gtransform_slice2( const GLU_complex *__restrict *__restrict gauge , struct site *__restrict lat , const GLU_complex *__restrict *__restrict gauge_up , const size_t t )
   @brief The gauge transform routine used in the coulomb gauge fixing routine
   @param gauge :: gauge transformation matrices for this timeslice
   @param lat :: lattice gauge field
   @param gauge_up :: gauge transformation matrices for the timeslice above
   @param t :: timeslice index
 **/
void
gtransform_slice2( const GLU_complex **gauge , 
		   struct site *lat , 
		   const GLU_complex **gauge_up ,
		   const size_t t ) ;


#endif
