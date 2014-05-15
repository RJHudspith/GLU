/*
    Copyright 2013 Renwick James Hudspith

    This file (pspace_landau.h) is part of GLU.

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
   @file pspace_landau.h
   @brief protype functions for the momentum space landau correction
 */
#ifndef GLU_PSPACE_LANDAU_H
#define GLU_PSPACE_LANDAU_H

/**
   @fn void correct_pspace_landau( struct site *__restrict A , const struct veclist *__restrict list , const int *__restrict in , const int DIMS )
   @brief performs the momentum space correction on our fields
   @param A :: Hermitian gauge fields in momentum space
   @param list :: the momentum lsit after the cutting procedure
   @param in :: the length of the momentum list , has dimension 1
   @param DIMS :: dimensions we are FFT-ing in

   Computes,
   \f[
   
      e^{ip_\mu/2} A_\mu(p)

   \f]
 */
void
correct_pspace_landau( struct site *__restrict A ,
		       const struct veclist *__restrict list ,
		       const int *__restrict in ,
		       const int DIMS ) ;

#endif
