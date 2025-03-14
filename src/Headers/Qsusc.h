/*
Copyright 2013-2025 Renwick James Hudspith

    This file (Qsusc.h) is part of GLU.

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
   @file Qsusc.h
   @brief prototype definitions for the Qsusc correlator
 */
#ifndef GLU_QSUSC_H
#define GLU_QSUSC_H

/**
   @fn int compute_Qsusc_step( struct site *__restrict lat , const struct cut_info CUTINFO , const struct sm_info SMINFO )
   @brief computes the topological chare correlation function
   @param lat :: lattice gauge fields
   @param CUTINFO :: uses the max_mom value for the maximum r^2 used
   @param SMINFO :: smearing information
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
compute_Qsusc_step( struct site *__restrict lat ,
		    const struct cut_info CUTINFO , 
		    const struct sm_info SMINFO ) ;

#endif
