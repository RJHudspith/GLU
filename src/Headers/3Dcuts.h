/*
Copyright 2013-2025 Renwick James Hudspith

    This file (3Dcuts.h) is part of GLU.

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
   @file 3Dcuts.h
   @brief Instantaneous temporal and spatial gluon propagator measurements
   @ingroup Cuts
 */
#ifndef GLU_CUTS_SPATIAL
#define GLU_CUTS_SPATIAL

/**
   @fn int cuts_spatial ( struct site *__restrict A ,
                          const struct cut_info CUTINFO ) ;

   @brief This code computes the instantaneous spatial and temporal gluon propagators this is for A - fields in Coulomb gauge.
   @param A :: The lattice field
   @param CUTINFO :: The general cutting information. i.e. Cut-type, Max and Min ... etc .

   As an added bonus it also shifts the result to reproduce the correct, 
   momentum-space coulomb condition. Being p.A == 0.

   @return #GLU_FAILURE or #GLU_SUCCESS
 */
int 
cuts_spatial ( struct site *__restrict A ,
	       const struct cut_info CUTINFO ) ;

#endif
