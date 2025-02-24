/*
Copyright 2013-2025 Renwick James Hudspith

    This file (POLY.h) is part of GLU.

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
   @file POLY.h
   @brief polyakov loop measurements, coulomb gauge fixed static potential meas
 */
#ifndef GLU_POLY_H
#define GLU_POLY_H

/**
   @fn int Coul_staticpot( struct site *lat , const struct cut_info CUTINFO , const struct sm_info SMINFO )
   @brief static potential calculator
   @param lat :: lattice gauge fields
   @param CUTINFO :: cutting information such as what type of cut to perform
   @param SMINFO :: do we want to smear this thing?
   @warning gauge configurations must be in Coulomb gauge, calls the smearing routines
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
Coul_staticpot( struct site *lat , 
		const struct cut_info CUTINFO ,
		const struct sm_info SMINFO ) ;

/**
   @fn double complex poly( const struct site *lat , size_t dir )
   @brief computes the polyakov loop in the direction "dir"
   @param lat :: lattice fields
   @param dir :: direction to measure in
   The measurement starts at the logical first sub-hypercube and traverses in
   the direction "dir" taking the product of the matrices from that point.

   @warning prints results to stdout. changes dir if not suitable
 **/
double complex 
poly( const struct site *lat , 
      size_t dir ) ;

/**
   @fn void poly_th( double *red , const struct site *lat , size_t dir ) ;
   @brief threaded polyakov loop evaluation
   @param red :: reduction array
   @param lat :: lattice gauge links
   @param dir :: direction to measure in between 0 and ND-1
 */
void
poly_th( double *red ,
	 const struct site *lat , 
	 size_t dir ) ;

#endif
