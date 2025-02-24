/*
Copyright 2013-2025 Renwick James Hudspith

    This file (gtrans.c) is part of GLU.

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
   @file gtrans.c
   @brief gauge transformations overwriting the lattice field
*/
#include "Mainfile.h"

#if (defined HAVE_IMMINTRIN_H) && !(defined SINGLE_PREC)

#include <immintrin.h>
#include "SSE2_OPS.h"

// inlined gauge transformation b = a.b.c^\dagger
static inline void
inline_gtransform_local( const __m128d *__restrict a ,
			 __m128d *__restrict b ,
			 const __m128d *__restrict c )
{
#if NC == 3
  register __m128d A = SSE2_MUL_CONJ( *(c+0) , *(b+0) ) ;
  register __m128d B = SSE2_MUL_CONJ( *(c+1) , *(b+1) ) ;
  register __m128d C = SSE2_MUL_CONJ( *(c+2) , *(b+2) ) ;
  register __m128d D = SSE2_MUL_CONJ( *(c+0) , *(b+3) ) ;
  register __m128d E = SSE2_MUL_CONJ( *(c+1) , *(b+4) ) ;
  register __m128d F = SSE2_MUL_CONJ( *(c+2) , *(b+5) ) ;

  A = _mm_add_pd( A , B ) ;
  D = _mm_add_pd( D , E ) ;
  const __m128d bp0 = _mm_add_pd( C , A ) ;
  const __m128d bp1 = _mm_add_pd( D , F ) ;

  A = SSE2_MUL_CONJ( *(c+0) , *(b+6) ) ;
  B = SSE2_MUL_CONJ( *(c+1) , *(b+7) ) ;
  D = SSE2_MUL_CONJ( *(c+3) , *(b+0) ) ;
  E = SSE2_MUL_CONJ( *(c+4) , *(b+1) ) ;
  C = SSE2_MUL_CONJ( *(c+2) , *(b+8) ) ;
  F = SSE2_MUL_CONJ( *(c+5) , *(b+2) ) ;

  A = _mm_add_pd( A , B ) ;
  D = _mm_add_pd( D , E ) ;
  const __m128d bp2 = _mm_add_pd( C , A ) ;
  const __m128d bp3 = _mm_add_pd( F , D ) ;

  A = SSE2_MUL_CONJ( *(c+3) , *(b+3) ) ;
  B = SSE2_MUL_CONJ( *(c+4) , *(b+4) ) ;
  D = SSE2_MUL_CONJ( *(c+3) , *(b+6) ) ;
  E = SSE2_MUL_CONJ( *(c+4) , *(b+7) ) ;
  C = SSE2_MUL_CONJ( *(c+5) , *(b+5) ) ;
  F = SSE2_MUL_CONJ( *(c+5) , *(b+8) ) ;

  A = _mm_add_pd( A , B ) ;
  D = _mm_add_pd( D , E ) ;
  const __m128d bp4 = _mm_add_pd( C , A ) ;
  const __m128d bp5 = _mm_add_pd( D , F ) ;

  A = SSE2_MULCONJ( bp0 , *( a + 0 ) ) ;
  B = SSE2_MULCONJ( bp1 , *( a + 1 ) ) ;
  D = SSE2_MULCONJ( bp3 , *( a + 0 ) ) ;
  E = SSE2_MULCONJ( bp4 , *( a + 1 ) ) ;
  C = SSE2_MULCONJ( bp2 , *( a + 2 ) ) ;
  F = SSE2_MULCONJ( bp5 , *( a + 2 ) ) ;

  A = _mm_add_pd( A , B ) ;
  D = _mm_add_pd( D , E ) ;
  *( b + 0 ) = _mm_add_pd( A , C ) ;
  *( b + 1 ) = _mm_add_pd( D , F ) ;

  A = SSE2_MULCONJ( bp0 , *( a + 3 ) ) ;
  B = SSE2_MULCONJ( bp1 , *( a + 4 ) ) ;
  D = SSE2_MULCONJ( bp3 , *( a + 3 ) ) ;
  E = SSE2_MULCONJ( bp4 , *( a + 4 ) ) ;
  C = SSE2_MULCONJ( bp2 , *( a + 5 ) ) ;
  F = SSE2_MULCONJ( bp5 , *( a + 5 ) ) ;
  A = _mm_add_pd( A , B ) ;
  D = _mm_add_pd( D , E ) ;
  *( b + 3 ) = _mm_add_pd( A , C ) ;
  *( b + 4 ) = _mm_add_pd( D , F ) ;

  A = SSE2_MULCONJ( bp0 , *( a + 6 ) ) ;
  B = SSE2_MULCONJ( bp1 , *( a + 7 ) ) ;
  D = SSE2_MULCONJ( bp3 , *( a + 6 ) ) ;
  E = SSE2_MULCONJ( bp4 , *( a + 7 ) ) ;
  C = SSE2_MULCONJ( bp2 , *( a + 8 ) ) ;
  F = SSE2_MULCONJ( bp5 , *( a + 8 ) ) ;
  A = _mm_add_pd( A , B ) ;
  D = _mm_add_pd( D , E ) ;
  *( b + 6 ) = _mm_add_pd( A , C ) ;
  *( b + 7 ) = _mm_add_pd( D , F ) ;

  A = SSE2_MUL( *( b + 3 ) , *( b + 7 ) ) ;
  B = SSE2_MUL( *( b + 4 ) , *( b + 6 ) ) ;
  C = SSE2_MUL( *( b + 1 ) , *( b + 6 ) ) ;
  D = SSE2_MUL( *( b + 0 ) , *( b + 7 ) ) ;
  E = SSE2_MUL( *( b + 0 ) , *( b + 4 ) ) ;
  F = SSE2_MUL( *( b + 1 ) , *( b + 3 ) ) ;

  A = _mm_sub_pd( A , B ) ;
  C = _mm_sub_pd( C , D ) ;
  E = _mm_sub_pd( E , F ) ;

  // complete
  *( b + 2 ) = SSE2_CONJ( A ) ;
  *( b + 5 ) = SSE2_CONJ( C ) ;
  *( b + 8 ) = SSE2_CONJ( E ) ;
#elif NC == 2
  // similar to above
  register const __m128d bp0 = _mm_add_pd( SSE2_MUL_CONJ( *(c+0) , *(b+0) ) ,
		    SSE2_MUL_CONJ( *(c+1) , *(b+1) ) ) ;
  register const __m128d bp1 = _mm_add_pd( SSE2_MUL_CONJ( *(c+0) , *(b+2) ) ,
		    SSE2_MUL_CONJ( *(c+1) , *(b+3) ) ) ;
  *(b+0) = _mm_add_pd( SSE2_MULCONJ( bp0 , *(a+0) ) ,
		       SSE2_MULCONJ( bp1 , *(a+1) ) ) ; 
  *(b+2) = _mm_add_pd( SSE2_MULCONJ( bp0 , *(a+2) ) ,
		       SSE2_MULCONJ( bp1 , *(a+3) ) ) ; 
  *(b+1) = SSE_FLIP( SSE2_CONJ( *( b + 2 ) ) ) ; 
  *(b+3) = SSE2_CONJ( *( b + 0 ) ) ;
#else
  // standard gauge transform
  GLU_complex temp[ NCNC ] GLUalign ;
  multab_dag_suNC( temp , (void*)b , (void*)c ) ;
  multab_suNC( (void*)b , (void*)a , temp ) ; 
#endif
  return ;
}
#else
// slower version
static inline void
inline_gtransform_local( const GLU_complex *__restrict a ,
			 GLU_complex *__restrict b ,
			 const GLU_complex *__restrict c )
{
  // standard gauge transform
  GLU_complex temp[ NCNC ] GLUalign ;
  multab_dag_suNC( temp , b , c ) ;
  multab_suNC( b , a , temp ) ; 
  return ;
}
#endif

// external version
void
gtransform_local( const GLU_complex *__restrict a ,
		  GLU_complex *__restrict b ,
		  const GLU_complex *__restrict c )
{
  inline_gtransform_local( (const void*)a , (void*)b , (const void*)c ) ;
  return ;
}

//gauge_transform lattice-wide
void 
gtransform( struct site *__restrict lat ,
	    const GLU_complex *__restrict *__restrict gauge )
{
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    #if ND == 4
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[0] , (const void*)gauge[lat[i].neighbor[0]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[1] , (const void*)gauge[lat[i].neighbor[1]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[2] , (const void*)gauge[lat[i].neighbor[2]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[3] , (const void*)gauge[lat[i].neighbor[3]] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[mu] , (const void*)gauge[lat[i].neighbor[mu]] ) ;
    }
    #endif
  } 
  return ;
}

//gauge_transform lattice-wide
void 
gtransform_th( struct site *__restrict lat ,
	       const GLU_complex *__restrict *__restrict gauge )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    #if ND == 4
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[0] , (const void*)gauge[lat[i].neighbor[0]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[1] , (const void*)gauge[lat[i].neighbor[1]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[2] , (const void*)gauge[lat[i].neighbor[2]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[3] , (const void*)gauge[lat[i].neighbor[3]] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      inline_gtransform_local( (const void*)gauge[i] , (void*)lat[i].O[mu] , (const void*)gauge[lat[i].neighbor[mu]] ) ;
    }
    #endif
  } 
  return ;
}

// gauge_transform for the Coulomb definition 
void
gtransform_slice( const GLU_complex *__restrict *__restrict gauge , 
		  struct site *__restrict lat , 
		  const GLU_complex *__restrict *__restrict gauge_up ,
		  const size_t t )
{
  size_t i ;
  const size_t slice = LCU  *  t ; 
#pragma omp parallel for private(i)
  for(  i = 0  ;  i < LCU  ;  i ++ ) {
    const size_t j = slice + i ;
    #if ND == 4   
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[0] , (const void*)gauge[lat[i].neighbor[0]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[1] , (const void*)gauge[lat[i].neighbor[1]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[2] , (const void*)gauge[lat[i].neighbor[2]] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND - 1  ; mu++ ) {
      inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[mu] , (const void*)gauge[lat[i].neighbor[mu]] ) ;
    }
    #endif
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[ND-1] , (const void*)gauge_up[i] ) ;
  }
  return ;
}

// gauge_transform for the Coulomb definition 
void
gtransform_slice_th( const GLU_complex *__restrict *__restrict gauge , 
		     struct site *__restrict lat , 
		     const GLU_complex *__restrict *__restrict gauge_up ,
		     const size_t t )
{
  size_t i ;
  const size_t slice = LCU  *  t ; 
#pragma omp for private(i)
  for(  i = 0  ;  i < LCU  ;  i ++ ) {
    const size_t j = slice + i ;
    #if ND == 4   
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[0] , (const void*)gauge[lat[i].neighbor[0]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[1] , (const void*)gauge[lat[i].neighbor[1]] ) ;
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[2] , (const void*)gauge[lat[i].neighbor[2]] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND - 1  ; mu++ ) {
      inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[mu] , (const void*)gauge[lat[i].neighbor[mu]] ) ;
    }
    #endif
    inline_gtransform_local( (const void*)gauge[i] , (void*)lat[j].O[ND-1] , (const void*)gauge_up[i] ) ;
  }
  return ;
}
