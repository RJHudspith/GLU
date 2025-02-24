/*
Copyright 2013-2025 Renwick James Hudspith

    This file (wflow.h) is part of GLU.

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
   @file wflow.h
   @brief prototype functions for the integration of gauge fields with a fictitious flow time
 */
#ifndef GLU_WFLOW_H
#define GLU_WFLOW_H

/**
   @fn int flow_RK3( struct site *__restrict lat , const size_t smiters , const int SIGN , const smearing_types SM_TYPE , const GLU_bool memcheap )
   @brief wilson flow using the Runge-Kutta used in Luescher's follow up paper and by BMW.
   @param lat :: lattice gauge field
   @param smiters :: number of smearing iterations
   @param SIGN :: of type #GLU_direction
   @param SM_TYPE :: will accept SM_STOUT or SM_LOG from #smearing_types
   @warning this one is the fastest but also the most memory expensive this is checked in the functions defined in GLU_memcheck.h
 **/
int
flow_RK3( struct site *__restrict lat , 
	  const size_t smiters ,
	  const int SIGN ,
	  const smearing_types SM_TYPE ,
	  const GLU_bool memcheap ) ;

#endif
