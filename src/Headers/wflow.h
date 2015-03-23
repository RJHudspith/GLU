/*
    Copyright 2013 Renwick James Hudspith

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
   @fn void flow4d_RK_fast( struct site *__restrict lat , const int smiters , const int DIR , const int SIGN , const int SM_TYPE )
   @brief wilson flow using the Runge-Kutta used in Luescher's follow up paper and by BMW.
   @param lat :: lattice gauge field
   @param lat :: lattice gauge fields 
   @param smiters :: number of smearing iterations
   @param DIR :: maximum direction of the smearing #GLU_smeardir
   @param SIGN :: of type #GLU_direction
   @param SM_TYPE :: will accept SM_APE , SM_STOUT or SM_LOG from #smearing_types
   @warning this one is the fastest but also the most memory expensive this is checked in the functions defined in GLU_memcheck.h
 **/
void 
flow4d_RK_fast( struct site *__restrict lat , 
		const int smiters ,
		const int DIR ,
		const int SIGN ,
		const int SM_TYPE ) ;

/**
   \fn void flow4d_RK_slow( struct site *__restrict lat , const int smiters , const int DIR , const int SIGN , const int SM_TYPE )
   \brief the terms slow and fast are really a matter of opinion, the slow one has fewer temporaries for the cost of (many) more memcpy's.
   @param lat :: lattice gauge field
   @param lat :: lattice gauge fields 
   @param smiters :: number of smearing iterations
   @param DIR :: maximum direction of the smearing #GLU_smeardir
   @param SIGN :: of type #GLU_direction
   @param SM_TYPE :: will accept SM_APE , SM_STOUT or SM_LOG from #smearing_types
   @warning selected depending on results in GLU_memcheck.h

   Uses the fact that we can iterate the smearing on a time-slice by time-slice basis just as we did in the 4D Hyp smearing codes (HYP_4D.h,4D_fast.h and New_HYP.h)
 **/
void 
flow4d_RK_slow( struct site *__restrict lat , 
		const int smiters ,
		const int DIR ,
		const int SIGN ,
		const int SM_TYPE ) ;

#endif
