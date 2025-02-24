/*
Copyright 2013-2025 Renwick James Hudspith

    This file (U1_obs.h) is part of GLU.

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
   @file U1_obs.h
   @brief U1 observable calculations
 */
#ifndef GLU_U1_OPS_H
#define GLU_U1_OPS_H

/**
   @fn void compute_U1_obs( const GLU_complex **U , const struct site *lat , const U1_meas meas )
   @brief this function is a wrapper for the U1 measurements
   @param U :: the U(1) non-compact field
   @param lat :: the lattice gauge field used for look up of neighbors
   @param meas :: the measurement being made of type #U1_meas
   it computes the non-compact and compact plaquettes by default. <br>
   if U1_RECTANGLE is defined, it outputs the compact rectangle <br>
   if U1_TOPOLOGICAL is defined it outputs the number of monopoles
   and the "dirac_sheet".
 **/
void
compute_U1_obs( const GLU_complex **U ,
		const struct site *lat ,
		const U1_meas meas ) ;

#endif
