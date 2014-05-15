/*
    Copyright 2013 Renwick James Hudspith

    This file (MOMggg.h) is part of GLU.

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
   @file MOMggg.h
   @brief gluon propagator and non-exceptional three point prototype function
   @ingroup Cuts
   non-exceptional scheme from <a href="http://arXiv.org/abs/1108.4806v1"> paper </a>
 */

#ifndef GLU_MOMGGG_h
#define GLU_MOMGGG_h

/**
   @fn int write_nonexceptional_g2g3( FILE *__restrict Ap , const struct site *__restrict A , const struct veclist *__restrict list , int num_mom[ 1 ] , const int nnmax )
   @brief computes the nonexceptional three point function and the gluon propagator
   @param Ap :: file we write out to
   @param A :: lie-algebra field
   @param list :: applicable momentum list after cuts 
   @param num_mom :: size of list
   @param nnmax :: maximum "p^2" allowed

   @return #GLU_SUCCESS or #GLU_FAILURE
 **/
int
write_nonexceptional_g2g3( FILE *__restrict Ap , 
			   const struct site *__restrict A , 
			   const struct veclist *__restrict list , 
			   int num_mom[ 1 ] , 
			   const int nnmax ) ;

#endif
