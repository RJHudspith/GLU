/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (GF_wrap.h) is part of GLU.

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
   @file GF_wrap.h
   @brief Funtion def for the wrapper of the gauge fixing code
 */
#ifndef GLU_GF_WRAP_H
#define GLU_GF_WRAP_H

/**
   @fn void GF_wrap( FILE *__restrict infile , struct site *__restrict lat , const struct gf_info GFINFO , const struct sm_info SMINFO , const struct head_data HEAD_DATA )
   @brief This is the wrapper that calls the gauge fixing functions.

   @param infile :: configuration file, be it HIREP or NERSC
   @param lat :: gauge field
   @param GFINFO :: general gauge fixing information from the input file
   @param SMINFO :: general smearing information from the input file
   @param HEAD_DATA :: information from the QCD header

   It includes the smeared-preconditioned and MAG-preconditioned
   types of gauge fixing.

   @return #GLU_FAILURE or #GLU_SUCCESS
 **/
int
GF_wrap( const char *__restrict infile , 
	 struct site *__restrict lat , 
	 const struct gf_info GFINFO , 
	 const struct sm_info SMINFO ,
	 const struct head_data HEAD_DATA ) ;

#endif
