/*
    Copyright 2013 Renwick James Hudspith

    This file (triplet_gen.h) is part of GLU.

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
   @file triplet_gen.h
   @brief prototype functions for generating triplets of external momenta for the non-exceptional scheme
   @ingroup Cuts
 */

#ifndef GLU_TRIPLET_GEN_H
#define GLU_TRIPLET_GEN_H

/**
   @fn int read_trip( int *__restrict trip , const size_t nnmax )
   @brief read the size of the triplets for each p^2
 */
int
read_trip( int *__restrict trip ,
	   const size_t nnmax ) ;

/**
   @fn int read_triplet_and_proj( int *__restrict *__restrict triplet , double *__restrict *__restrict proj , const int *__restrict *__restrict momentum , const size_t nnmax , const size_t count , const size_t nmom )
   @brief computes the triplet and projector
   @param triplet :: is of the form triplet[ triplet_index ][ 0,1,2 ] as there are three of them
   @param proj :: projector for obtaining the non-exceptional scalar three point function
   @param momentum :: a list of the momentum in integer rep of the form momentum[ index ][ 0,1,2 .. ND ]
   @param nnmax :: maximum \f$p^2\f$
   @param count :: the number of all triplets
   @param nmom :: the number of momentum after cutting
   @warning nnmax must be even
   @bug nnmax which is not even breaks this for some reason

   <br>
   Either reads/generates a file if "NOT_CONDOR_MODE" is defined, then computes
   both the triplet's places in the momentum configuration and the projector
   that we are using.
   <br>

   @return #GLU_SUCCESS or #GLU_FAILURE
 **/
int
read_triplet_and_proj( int *__restrict *__restrict triplet , 
		       double *__restrict *__restrict proj , 
		       const int *__restrict *__restrict momentum , 
		       const size_t nnmax , 
		       const size_t count ,
		       const size_t nmom ) ;

#endif
