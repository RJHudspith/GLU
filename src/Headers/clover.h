/*
Copyright 2013-2025 Renwick James Hudspith

    This file (clover.h) is part of GLU.

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
   @file clover.h 
   @brief Calculation of the naive field strength "G_{\mu\nu}" and topological charge
 */
#ifndef GLU_CLOVER_H
#define GLU_CLOVER_H

/**
   @fn void compute_Gmunu_th( double *red , const struct site *lat )
   @brief Gmunu calculation to be called within a parallel environment
   @param red :: reduction array
   @param lat :: lattice gauge links
 */
void
compute_Gmunu_th( double *red , 
		  const struct site *lat ) ;

/**
   @fn void compute_Gmunu_array( GLU_complex *__restrict *__restrict qtop , const struct site *__restrict lat )
   @brief \f$ O(a^4) \f$ tree-improved field strength tensor from <a href="http://arxiv.org/abs/hep-lat/0203008"> paper </a>
   @param qtop :: Naive topological charge matrix
   @param lat :: lattice field
 */
void
compute_Gmunu_array( GLU_complex *__restrict qtop ,
		     const struct site *__restrict lat ) ;

#endif
