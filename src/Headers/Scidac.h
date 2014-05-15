/*
    Copyright 2013 Renwick James Hudspith

    This file (Scidac.h) is part of GLU.

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
   @file Scidac.h
   @brief function prototypes for the reading of SCIDAC files
 */

#ifndef GLU_SCIDAC_H
#define GLU_SCIDAC_H

/**
   @fn int get_header_data_SCIDAC( FILE *infile , struct head_data *HEAD_DATA )
   @brief gets the configuration data from the Scidac or ILDG header
   @param infile :: configuration file
   @param HEAD_DATA :: contains the necessary header data

   @warning overwrites the struct Latt's dimensions
*/
int
get_header_data_SCIDAC( FILE *infile ,
			struct head_data *HEAD_DATA ) ;

/**
   @fn void write_header_ILDG( FILE *__restrict out )
   @brief function for writing out an ILDG header
   @param out :: output file
 */
void
write_header_ILDG( FILE *__restrict out ) ;

/**
   @fn void write_header_SCIDAC( FILE *__restrict out ) ;
   @brief function for writing out a SCIDAC configuration
   @param out :: the output file
 */
void
write_header_SCIDAC( FILE *__restrict out ) ;

/**
   @fn void write_trailing_header_SCIDAC( FILE *__restrict out , const uint32_t cksuma , const uint32_t cksumb ) 
   @brief writes the checksums for the configuration file
   @param cksuma :: the computed CRC checksum (rank29)
   @param cksumb :: the CRC checksum (rank31)
 */
void
write_trailing_header_SCIDAC( FILE *__restrict out ,
			      const uint32_t cksuma ,
			      const uint32_t cksumb ) ;

#endif
