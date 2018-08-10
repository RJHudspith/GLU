/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (Landau.h) is part of GLU.

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
   @file Landau.h
   @brief prototype functions for the Landau gauge fixing
 */
#ifndef GLU_LANDAU_H
#define GLU_LANDAU_H

/**
   @fn int grab_file( struct site *__restrict lat , GLU_complex *__restrict *__restrict gauge , FILE *__restrict infile )
   @brief reread the configuration file and perform a random transform using gauge
   @param lat :: lattice links, overwritten
   @param gauge :: gauge transformation matrices
   @param infile :: input file
   
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
grab_file( struct site *lat , 
	   GLU_complex **gauge , 
	   const char *infile ) ;

/**
   @fn size_t Landau( struct site *lat , const double accuracy , const size_t iter , FILE *infile )
   @brief fast landau gauge fixing
   @param lat :: gauge field
   @param accuracy :: accuracy we wish to converge our algorithm to
   @param iter :: maximum number of iterations before we restart
   @param infile :: need the configuration file in case we have to start over
   our fastest smearing routine, memory expensive
   @return the number of iterations the routine took
 **/
size_t
Landau( struct site *lat ,
	const double accuracy ,
	const size_t iter ,
	const char *infile ) ;

#endif
