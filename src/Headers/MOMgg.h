/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (MOMgg.h) is part of GLU.

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
   @file MOMgg.h
   @brief calculation of the gluon propagator and the exceptional three point function 
   @ingroup Cuts
   The projection is taken from <a href="http://arXiv:hep-ph/abs/0007088v1"> paper </a>
*/

#ifndef GLU_MOMGG_H
#define GLU_MOMGG_H

/**
   @fn int write_exceptional_g2g3_MOMgg( FILE *__restrict Ap , const struct site *__restrict A , const struct veclist *__restrict list , int num_mom[1] )
   @brief gluon propagator and exceptional three point   
   @param Ap :: file we write out to
   @param A :: lie-algebra field
   @param list :: applicable momentum list after cuts 
   @param num_mom :: maximum number of momenta 
   Includes the calculation of the gluonic two point correlation
   functions and the three point function with the MOMgg kinematics/projector.
   
   @return #GLU_SUCCESS or #GLU_FAILURE
 **/
int
write_exceptional_g2g3_MOMgg( FILE *__restrict Ap , 
			      const struct site *__restrict A , 
			      const struct veclist *__restrict list , 
			      int num_mom[1] ) ;

#endif
