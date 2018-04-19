/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (crc.h) is part of GLU.

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
   @file crc.h
   @brief the crc for the scidac checksum
 */
#ifndef GLU_CRC_H
#define GLU_CRC_H

/**
   @fn uint32_t crc32( uint32_t crc , const unsigned char *buf , size_t len )
   @brief takes the crc checksum of buf 
   @param crc :: previous CRC checksum
   @param buf :: the data having its CRC taken
   @param len :: number of bytes
 */
uint32_t
crc32( uint32_t crc , 
       const unsigned char *buf , 
       size_t len ) ;

/**
   @fn void CKSUM_ADD( void *memptr , const uint32_t nbytes )
   @brief BQCD's crc accumulator
   @param memptr :: start of the memory that we crc
   @param nbytes :: total number of bytes allocated in memptr
 */
void 
CKSUM_ADD( void *memptr , 
	   const uint32_t nbytes ) ;

/**
   @fn void CKSUM_GET( uint32_t *total_crc, uint32_t *total_bytes ) ;
   @brief BQCD's crc
 */
void 
CKSUM_GET( uint32_t *total_crc, 
	   uint32_t *total_bytes ) ;

#endif
