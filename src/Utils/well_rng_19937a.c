/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (well_rng_19937a.c) is part of GLU.

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
   @file well_rng_19937a.c
   @brief L'Ecuyer,Panneton and Matsumoto's WELL random number generator
 */

/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */

#include <stdint.h>
#include <stdlib.h>

#define W 32
#define R 624
#define P 31
#define MASKU (0xffffffffU>>(W-P))
#define MASKL (~MASKU)
#define M1 70
#define M2 179
#define M3 449

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT1(v) v
#define MAT3POS(t,v) (v>>t)

/* To obtain the WELL19937c, uncomment the following line */
/* 
#define TEMPERING                                      
static uint32_t y;
#define TEMPERB 0xe46e1700U
#define TEMPERC 0x9b868000U
*/

static uint32_t *table , z0  , z1 , z2 , state_i = 0 ;

#define V0            table[state_i]
#define VM1Over       table[state_i+M1-R]
#define VM1           table[state_i+M1]
#define VM2Over       table[state_i+M2-R]
#define VM2           table[state_i+M2]
#define VM3Over       table[state_i+M3-R]
#define VM3           table[state_i+M3]
#define VRm1          table[state_i-1]
#define VRm1Under     table[state_i+R-1]
#define VRm2          table[state_i-2]
#define VRm2Under     table[state_i+R-2]

#define newV0         table[state_i-1]
#define newV0Under    table[state_i-1+R]
#define newV1         table[state_i]
#define newVRm1       table[state_i-2]
#define newVRm1Under  table[state_i-2+R]

#define FACT 2.32830643653869628906e-10

static double case_1 (void);
static double case_2 (void);
static double case_3 (void);
static double case_4 (void);
static double case_5 (void);
static double case_6 (void);
double (*WELLRNG19937a) (void);

// set the table
void
GLU_set_WELL19937_table( const uint32_t seed )
{
  table = ( uint32_t* )malloc( 624 * sizeof( uint32_t ) ) ;
  table[0] = seed ;
  int i ;
  for( i = 1 ; i < 624 ; i++ ) {
    table[i] = ( 1812433253UL * ( table[i-1] ^ ( table[i-1] >> 30)) + i ) ;
  }
  // set the state and point at case_1
  state_i = 0 ;
  WELLRNG19937a = case_1 ;
  return ;
}

// free the table
void
GLU_free_WELL19937_table( void )
{
  free( table ) ;
}

static double case_1 (void)
{
   // state_i == 0
   z0 = (VRm1Under & MASKL) | (VRm2Under & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3);
   newV1      = z1 ^ z2;
   newV0Under = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i = R - 1;
   WELLRNG19937a = case_3;
#ifdef TEMPERING
   y = table[state_i] ^ ((table[state_i] << 7) & TEMPERB);
   y =              y ^ ((             y << 15) & TEMPERC);
   return ((double) y * FACT);
#else
   return ((double) table[state_i] * FACT);
#endif
}

static double case_2 (void)
{
   // state_i == 1
   z0 = (VRm1 & MASKL) | (VRm2Under & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i = 0;
   WELLRNG19937a = case_1;
#ifdef TEMPERING
   y = table[state_i] ^ ((table[state_i] << 7) & TEMPERB);
   y =              y ^ ((             y << 15) & TEMPERC);
   return ((double) y * FACT);
#else
   return ((double) table[state_i] * FACT);
#endif
}

static double case_3 (void)
{
   // state_i+M1 >= R
   z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1Over);
   z2 = MAT3POS (9, VM2Over) ^ MAT0POS (1, VM3Over);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i--;
   if (state_i + M1 < R)
      WELLRNG19937a = case_5;
#ifdef TEMPERING
   y = table[state_i] ^ ((table[state_i] << 7) & TEMPERB);
   y =              y ^ ((             y << 15) & TEMPERC);
   return ((double) y * FACT);
#else
   return ((double) table[state_i] * FACT);
#endif
}

static double case_4 (void)
{
   // state_i+M3 >= R
   z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3Over);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i--;
   if (state_i + M3 < R)
      WELLRNG19937a = case_6;
#ifdef TEMPERING
   y = table[state_i] ^ ((table[state_i] << 7) & TEMPERB);
   y =              y ^ ((             y << 15) & TEMPERC);
   return ((double) y * FACT);
#else
   return ((double) table[state_i] * FACT);
#endif
}

static double case_5 (void)
{
   // state_i+M2 >= R
   z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2Over) ^ MAT0POS (1, VM3Over);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i--;
   if (state_i + M2 < R)
      WELLRNG19937a = case_4;
#ifdef TEMPERING
   y = table[state_i] ^ ((table[state_i] << 7) & TEMPERB);
   y =              y ^ ((             y << 15) & TEMPERC);
   return ((double) y * FACT);
#else
   return ((double) table[state_i] * FACT);
#endif
}

static double case_6 (void)
{
   // 2 <= state_i <= (R - M3 - 1)
   z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i--;
   if (state_i == 1)
      WELLRNG19937a = case_2;
#ifdef TEMPERING
   y = table[state_i] ^ ((table[state_i] << 7) & TEMPERB);
   y =              y ^ ((             y << 15) & TEMPERC);
   return ((double) y * FACT);
#else
   return ((double) table[state_i] * FACT);
#endif
}

// undefs
#undef W
#undef R
#undef P
#undef MASKU
#undef MASKL
#undef M1
#undef M2
#undef M3
#undef MAT0POS
#undef MAT0NEG
#undef MAT1
#undef MAT3POS
#undef V0
#undef VM1Over
#undef VM1     
#undef VM2Over
#undef VM2
#undef VM3Over
#undef VM3
#undef VRm1
#undef VRm1Under
#undef VRm2
#undef VRm2Under
#undef newV0
#undef newV0Under
#undef newV1
#undef newVRm1
#undef newVRm1Under
#ifdef TEMPERING
 #undef TEMPERING
 #undef TEMPERB
 #undef TEMPERC
#endif
