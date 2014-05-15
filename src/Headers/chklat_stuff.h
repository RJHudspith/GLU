/*
    Copyright 2013 Renwick James Hudspith

    This file (chklat_stuff.h) is part of GLU.

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
   @file chklat_stuff.h
   @brief prototype functions from chklat.c that I found useful <a href="linkURL">http://lattices.qcdoc.bnl.gov/formatDescription/chklat_32.c </a>

   mostly for getting info from the header
 */

#ifndef GLU_CHKLAT_STUFF
#define GLU_CHKLAT_STUFF

/**
   @fn struct QCDheader * get_header( FILE *__restrict in )
   @brief Chklat's QCDheader reader

   @param in :: NERSC format file being read
 **/
struct QCDheader * 
get_header( FILE *__restrict in ) ;

/**
   @fn int get_string( char *s , struct QCDheader *hdr , char **q )
   @brief Gets a string from the header 

   @param s :: The string tag we are looking for
   @param hdr :: the NERSC header we have read in
   @param q :: the string being returned
 **/
int 
get_string( char *s , 
	    struct QCDheader *hdr ,
	    char **q ) ;

/**
   @fn int get_uint32_t( char *s , struct QCDheader *hdr , uint32_t *q )
   @brief Gets a uint32_t from the header 

   @param s :: The string tag we are looking for
   @param hdr :: the NERSC header we have read in
   @param q :: the int being returned

   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int 
get_uint32_t( char *s , 
	      struct QCDheader *hdr ,
	      uint32_t *q ) ;

/**
   @fn get_int( char *s , struct QCDheader *hdr , int *q )
   @brief Gets an int from the header 

   @param s :: The string tag we are looking for
   @param hdr :: the NERSC header we have read in
   @param q :: the int being returned

   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int 
get_int( char *s , 
	 struct QCDheader *hdr ,
	 int *q ) ;

/**
   @fn get_float( char *s , struct QCDheader *hdr , float *q ) ;
   @brief Gets a float from the header 

   @param s :: The string tag we are looking for
   @param hdr :: the NERSC header we have read in
   @param q :: the int being returned

   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int 
get_float( char *s , 
	   struct QCDheader *hdr ,
	   float *q ) ;

/**
   @fn int skip_hdr( FILE *__restrict file ) ;
   @brief Skips the header completely, dangeous
   @param file :: NERSC configuration file
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int 
skip_hdr( FILE *__restrict file ) ;

#endif
