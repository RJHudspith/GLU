/*
Copyright 2013-2025 Renwick James Hudspith

    This file (cuts.h) is part of GLU.

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
   @file cuts.h
   @brief prototype function def for the cutting procedures
   @ingroup Cuts
 */
#ifndef GLU_CUTS_H
#define GLU_CUTS_H

/**
   @fn int cuts_struct( struct site *__restrict A , const struct cut_info CUTINFO )
   @brief cutting procedure generating our Fourier-transformed A(p) lie fields.
   @param A :: Our Fourier transformed gauge fields
   @param CUTINFO :: General cutting information

   Same as the old routine but with all of the cutting information
   squirreled away in CUTINFO..

   @warning overwrites the stored gauge links in A

   @return #GLU_SUCCESS or #GLU_FAILURE
 **/
int 
cuts_struct( struct site *__restrict A ,
	     const struct cut_info CUTINFO ) ;

#endif
