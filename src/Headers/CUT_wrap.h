/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (CUT_wrap.h) is part of GLU.

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
   @file CUT_wrap.h
   @brief Wrapper for the cutting procedure
 */

#ifndef GLU_CUT_WRAP
#define GLU_CUT_WRAP

/**
   @fn void cuts_wrap_struct( struct site *__restrict lat , const struct cut_info CUTINFO , const struct sm_info SMINFO )
   @brief Wrapper function

   @param lat :: Gauge fields
   @param CUTINFO :: cutting information
   @param SMINFO :: smearing information
 */
void
cuts_wrap_struct( struct site *__restrict lat , 
		  const struct cut_info CUTINFO ,
		  const struct sm_info SMINFO ) ;

#endif
