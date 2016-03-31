/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (adaptive_flow.h) is part of GLU.

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
   @file adaptive_flow.h
   @brief prototype function defs for the adaptive routine
 */

#ifndef GLU_ADAPTIVE_FLOW_H
#define GLU_ADAPTIVE_FLOW_H

/**
   @fn int flow4d_adaptive_RK( struct site *__restrict lat , const int smiters , const int DIR , const int SIGN , const smearing_types SM_TYPE ) ;
   @brief this code performs a two-step adaptive RK integration of the flow equation
   @param lat :: lattice gauge field
   @param lat :: lattice gauge fields 
   @param smiters :: number of smearing iterations
   @param DIR :: maximum direction of the smearing #GLU_smeardir
   @param SIGN :: of type #GLU_direction
   @param SM_TYPE :: will accept SM_APE , SM_STOUT or SM_LOG from #smearing_types
   <br>
   it is expensive in memory (two lattice fields needed to compare) and uses the same static functions as void flow4d_RK_fast( struct site *lat , const int smiters , const int DIR , const int SIGN , const smearing_types SM_TYPE ) and so is comparable in speed. The code also shortens the time step when we get close to WFLOW_STOP, so we can measure W0 accurately.
 **/
int 
flow4d_adaptive_RK( struct site *__restrict lat , 
		    const size_t smiters ,
		    const size_t DIR ,
		    const int SIGN ,
		    const smearing_types SM_TYPE ) ;


#endif
