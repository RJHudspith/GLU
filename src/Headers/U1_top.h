/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (U1_top.h) is part of GLU.

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
   @file U1_top.h
   @brief U1 topological measures prototype function
 */
#ifndef GLU_U1_TOP_H
#define GLU_U1_TOP_H

/**
   @fn void U1_topological( int *__restrict monopole , int *__restrict d_sheet , double *__restrict qtop , const GLU_real *__restrict *__restrict O  )
   @brief measures some of the more elementary topological defects from the noncompact U(1) gauge field
   @param monopole :: number of monopoles
   @param d_sheet :: the "dirac sheet", the sum total of the non-zero windings around the elementary plaquette
   @param qtop :: the simplest measure of the U(1) naive gauge topological charge
   @param O :: the noncompact U(1) gauge field 
 */
void 
U1_topological( int *__restrict monopole , 
		int *__restrict d_sheet , 
		double *__restrict qtop ,
		const GLU_real *__restrict *__restrict O ) ;

#endif
