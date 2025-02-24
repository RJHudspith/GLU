/*
Copyright 2013-2025 Renwick James Hudspith

    This file (glueprop.h) is part of GLU.

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
   @file glueprop.h
   @brief prototype function def for gluon propagator measurements
   @ingroup Cuts
 */
#ifndef GLU_GLUEPROP_H
#define GLU_GLUEPROP_H

/**
   @fn int compute_gluon_prop( FILE *__restrict Ap , const struct site *__restrict A , const struct veclist *__restrict list , size_t num_mom[1] )
   @brief computes the transverse and longitudinal gluon propagator scalar functions
   @param Ap :: file being written out to
   @param A :: Momentum space gluon fields
   @param list :: momentum list after cutting
   @param num_mom :: length of the momentum list

   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
compute_gluon_prop( FILE *__restrict Ap , 
		    const struct site *__restrict A ,
		    const struct veclist *__restrict list ,
		    size_t num_mom[1] ) ;

#endif
