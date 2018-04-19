/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (write_headers.h) is part of GLU.

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
   @file write_headers.h
   @brief prototypes for writing the header information of MILC and NERSC type configuration files
 */
#ifndef GLU_WRITE_HEADERS_H
#define GLU_WRITE_HEADERS_H

/**
   @fn void write_header_MILC( FILE *__restrict out , const uint32_t milc_cksum29 , const uint32_t milc_cksum31 ) ;
   @brief writes the MILC configuration header to the start of out
   @param out :: the output file being written to
   @param milc_cksum29 :: the bit-shifted by index (mod 29) checksum
   @param milc_cksum29 :: the bit-shifted by index (mod 29) checksum
 */
void
write_header_MILC( FILE *__restrict out ,
		   const uint32_t milc_cksum29 ,
		   const uint32_t milc_cksum31 ) ;

/**
   @fn void write_header_NERSC( FILE *__restrict out , const double tr , const double plaq , const uint32_t chksum , const char *details , const int type ) ;
   @brief writes the NERSC header format to the top of "out"
   @param out :: file to be written in
   @param tr :: average trace of the lattice links
   @param plaq :: average plaquette of the lattice fields
   @param chksum :: NERSC checksum is a sum of the data cast to uint32_t
   @param details :: what goes in the INFO part of the header
   @param type :: can be #OUTPUT_NCxNC,#OUTPUT_GAUGE or #OUTPUT_SMALL
 */
void
write_header_NERSC( FILE *__restrict out ,
		    const double tr , 
		    const double plaq ,
		    const uint32_t chksum ,
		    const char *details ,
		    const int type ) ;

#endif
