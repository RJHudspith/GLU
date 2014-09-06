/*
    Copyright 2013 Renwick James Hudspith

    This file (GLU_bswap.c) is part of GLU.

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
   @file GLU_bswap.c
   @brief byte swapping utilities and routines
   
   All of the external routines require an array to be swapped
 */

#include <stdint.h>

// These three functions are from LLVM, credit to them, I didn't know
// gcc had an intrinsic for this
static inline uint16_t 
SwapByteOrder_16( uint16_t value ) {
#if defined(_MSC_VER) && !defined(_DEBUG)
  // The DLL version of the runtime lacks these functions (bug!?), but in a
  // release build they're replaced with BSWAP instructions anyway.
  return _byteswap_ushort(value);
#else
  uint16_t Hi = value << 8;
  uint16_t Lo = value >> 8;
  return Hi | Lo;
#endif
}

// SwapByteOrder_32 - This function returns a byte-swapped representation of
// the 32-bit argument.
static inline uint32_t 
SwapByteOrder_32( uint32_t value ) {
#if defined(__llvm__) || (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)) && !defined(__ICC)
  return __builtin_bswap32(value);
#elif defined(_MSC_VER) && !defined(_DEBUG)
  return _byteswap_ulong(value);
#else
  uint32_t Byte0 = value & 0x000000FF ;
  uint32_t Byte1 = value & 0x0000FF00 ;
  uint32_t Byte2 = value & 0x00FF0000 ;
  uint32_t Byte3 = value & 0xFF000000 ;
  return (Byte0 << 24) | (Byte1 << 8) | (Byte2 >> 8) | (Byte3 >> 24) ;
#endif
}

// SwapByteOrder_64 - This function returns a byte-swapped representation of
// the 64-bit argument.
static inline uint64_t 
SwapByteOrder_64( uint64_t value ) {
#if defined(__llvm__) || (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)) && !defined(__ICC)
   return __builtin_bswap64(value) ;
#elif defined(_MSC_VER) && !defined(_DEBUG)
   return _byteswap_uint64(value) ;
#else
   uint64_t Hi = SwapByteOrder_32( (uint32_t)( value ) ) ;
   uint32_t Lo = SwapByteOrder_32( (uint32_t)( value >> 32 ) ) ;
   return (Hi << 32) | Lo ;
#endif
}

// this is for 16 bits, i.e. shorts
void
bswap_16( const int n , void *u ) 
{
  uint16_t *T = u ;
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < n ; i++ ) {
    *( T + i ) = SwapByteOrder_16( *( T + i ) ) ;
  }
  return ;
} 

// this is for 32 bits, i.e. float and int
void
bswap_32( const int n , void *u ) 
{
  uint32_t *T = u ;
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < n ; i++ ) {
    *( T + i ) = SwapByteOrder_32( *( T + i ) ) ;
  }
  return ;
} 

// this one is for 64 bits, double and double complex
void
bswap_64( const int n , void *u ) 
{
  uint64_t *T = u ;
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < n ; i++ ) {
    *( T + i ) = SwapByteOrder_64( *( T + i ) ) ;
  }
  return ;
}
