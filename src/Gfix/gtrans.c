/*
    Copyright 2013-2016 Renwick James Hudspith

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

static inline void
inline_gtransform_local( const GLU_complex *__restrict a ,
			 GLU_complex *__restrict b ,
			 const GLU_complex *__restrict c )
{
#if NC == 3
  // standard gauge transform
  __m128d *B = (__m128d*)b ;
  const __m128d *A = (const __m128d*)a ;
  const __m128d *C = (const __m128d*)c ;

  register __m128d c0 = *C ; C++ ; 
  register __m128d c1 = *C ; C++ ; 
  register __m128d c2 = *C ; C++ ;
  register __m128d c3 = *C ; C++ ; 
  register __m128d c4 = *C ; C++ ; 
  register __m128d c5 = *C ; // falls out of cache here!!

  register const __m128d bp0 = _mm_add_pd( SSE2_MUL_CONJ( c0 , *B ) ,
		    _mm_add_pd( SSE2_MUL_CONJ( c1 , *( B + 1 ) ) ,
				SSE2_MUL_CONJ( c2 , *( B + 2 ) ) ) ) ;
  register const __m128d bp1 = _mm_add_pd( SSE2_MUL_CONJ( c0 , *( B + 3 ) ) ,
		    _mm_add_pd( SSE2_MUL_CONJ( c1 , *( B + 4 ) ) ,
				SSE2_MUL_CONJ( c2 , *( B + 5 ) ) ) ) ;
  register const __m128d bp2 = _mm_add_pd( SSE2_MUL_CONJ( c0 , *( B + 6 ) ) ,
		    _mm_add_pd( SSE2_MUL_CONJ( c1 , *( B + 7 ) ) ,
				SSE2_MUL_CONJ( c2 , *( B + 8 ) ) ) ) ;
  // compute the second row
  register const __m128d bp3 = _mm_add_pd( SSE2_MUL_CONJ( c3 , *B ) ,
		    _mm_add_pd( SSE2_MUL_CONJ( c4 , *( B + 1 ) ) ,
				SSE2_MUL_CONJ( c5 , *( B + 2 ) ) ) ) ;
  register const __m128d bp4 = _mm_add_pd( SSE2_MUL_CONJ( c3 , *( B + 3 ) ) ,
		    _mm_add_pd( SSE2_MUL_CONJ( c4 , *( B + 4 ) ) ,
				SSE2_MUL_CONJ( c5 , *( B + 5 ) ) ) ) ;
  register const __m128d bp5 = _mm_add_pd( SSE2_MUL_CONJ( c3 , *( B + 6 ) ) ,
		    _mm_add_pd( SSE2_MUL_CONJ( c4 , *( B + 7 ) ) ,
				SSE2_MUL_CONJ( c5 , *( B + 8 ) ) ) ) ;
  // set the first two columns of B
  *( B + 0 ) = _mm_add_pd( SSE2_MULCONJ( bp0 , *( A + 0 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( bp1 , *( A + 1 ) ) ,
				       SSE2_MULCONJ( bp2 , *( A + 2 ) ) ) ) ;
  *( B + 1 ) = _mm_add_pd( SSE2_MULCONJ( bp3 , *( A + 0 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( bp4 , *( A + 1 ) ) ,
				       SSE2_MULCONJ( bp5 , *( A + 2 ) ) ) ) ;
  *( B + 3 ) = _mm_add_pd( SSE2_MULCONJ( bp0 , *( A + 3 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( bp1 , *( A + 4 ) ) ,
				       SSE2_MULCONJ( bp2 , *( A + 5 ) ) ) ) ;
  *( B + 4 ) = _mm_add_pd( SSE2_MULCONJ( bp3 , *( A + 3 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( bp4 , *( A + 4 ) ) ,
				       SSE2_MULCONJ( bp5 , *( A + 5 ) ) ) ) ;
  *( B + 6 ) = _mm_add_pd( SSE2_MULCONJ( bp0 , *( A + 6 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( bp1 , *( A + 7 ) ) ,
				       SSE2_MULCONJ( bp2 , *( A + 8 ) ) ) ) ;
  *( B + 7 ) = _mm_add_pd( SSE2_MULCONJ( bp3 , *( A + 6 ) ) ,
			   _mm_add_pd( SSE2_MULCONJ( bp4 , *( A + 7 ) ) ,
				       SSE2_MULCONJ( bp5 , *( A + 8 ) ) ) ) ;
  // complete
  *( B + 2 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( B + 3 ) , *( B + 7 ) ) ,
				      SSE2_MUL( *( B + 4 ) , *( B + 6 ) ) ) ) ;
  *( B + 5 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( B + 1 ) , *( B + 6 ) ) ,
				      SSE2_MUL( *( B + 0 ) , *( B + 7 ) ) ) ) ;
  *( B + 8 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( B + 0 ) , *( B + 4 ) ) ,
				      SSE2_MUL( *( B + 1 ) , *( B + 3 ) ) ) ) ;
#elif NC == 2
  __m128d *B = (__m128d*)b ;
  const __m128d *A = (const __m128d*)a ;
  const __m128d *C = (const __m128d*)c ;
  // similar to above
  register const __m128d bp0 = _mm_add_pd( SSE2_MUL_CONJ( *(C+0) , *(B+0) ) ,
		    SSE2_MUL_CONJ( *(C+1) , *(B+1) ) ) ;
  register const __m128d bp1 = _mm_add_pd( SSE2_MUL_CONJ( *(C+0) , *(B+2) ) ,
		    SSE2_MUL_CONJ( *(C+1) , *(B+3) ) ) ;
  *(B+0) = _mm_add_pd( SSE2_MULCONJ( bp0 , *(A+0) ) ,
		       SSE2_MULCONJ( bp1 , *(A+1) ) ) ; 
  *(B+2) = _mm_add_pd( SSE2_MULCONJ( bp0 , *(A+2) ) ,
		       SSE2_MULCONJ( bp1 , *(A+3) ) ) ; 
  *(B+1) = SSE_FLIP( SSE2_CONJ( *( B + 2 ) ) ) ; 
  *(B+3) = SSE2_CONJ( *( B + 0 ) ) ;
#else
  // standard gauge transform
  GLU_complex temp[ NCNC ] GLUalign ;
  multab_dag_suNC( temp , b , c ) ;
  multab_suNC( b , a , temp ) ; 
#endif
  return ;
}
#else
// slow version
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
  inline_gtransform_local( a , b , c ) ;
  return ;
}

//gauge_transform lattice-wide
void 
gtransform( struct site *__restrict lat ,
	    const GLU_complex *__restrict *__restrict gauge )
{
  size_t i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    #if ND == 4
    inline_gtransform_local( gauge[i] , lat[i].O[0] , gauge[lat[i].neighbor[0]] ) ;
    inline_gtransform_local( gauge[i] , lat[i].O[1] , gauge[lat[i].neighbor[1]] ) ;
    inline_gtransform_local( gauge[i] , lat[i].O[2] , gauge[lat[i].neighbor[2]] ) ;
    inline_gtransform_local( gauge[i] , lat[i].O[3] , gauge[lat[i].neighbor[3]] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      inline_gtransform_local( gauge[i] , lat[i].O[mu] , gauge[lat[i].neighbor[mu]] ) ;
    }
    #endif
  } 
  return ;
}

//gauge_transform lattice-wide
void 
gtransform2( struct site *__restrict lat ,
	     const GLU_complex *__restrict *__restrict gauge )
{
  size_t i ;
#pragma omp for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    #if ND == 4
    inline_gtransform_local( gauge[i] , lat[i].O[0] , gauge[lat[i].neighbor[0]] ) ;
    inline_gtransform_local( gauge[i] , lat[i].O[1] , gauge[lat[i].neighbor[1]] ) ;
    inline_gtransform_local( gauge[i] , lat[i].O[2] , gauge[lat[i].neighbor[2]] ) ;
    inline_gtransform_local( gauge[i] , lat[i].O[3] , gauge[lat[i].neighbor[3]] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      inline_gtransform_local( gauge[i] , lat[i].O[mu] , gauge[lat[i].neighbor[mu]] ) ;
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
  PFOR(  i = 0  ;  i < LCU  ;  i ++ ) {
    const size_t j = slice + i ;
    #if ND == 4   
    inline_gtransform_local( gauge[i] , lat[j].O[0] , gauge[lat[i].neighbor[0]] ) ;
    inline_gtransform_local( gauge[i] , lat[j].O[1] , gauge[lat[i].neighbor[1]] ) ;
    inline_gtransform_local( gauge[i] , lat[j].O[2] , gauge[lat[i].neighbor[2]] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND - 1  ; mu++ ) {
      inline_gtransform_local( gauge[i] , lat[j].O[mu] , gauge[lat[i].neighbor[mu]] ) ;
    }
    #endif
    inline_gtransform_local( gauge[i] , lat[j].O[ND-1] , gauge_up[i] ) ;
  }
  return ;
}

// gauge_transform for the Coulomb definition 
void
gtransform_slice2( const GLU_complex *__restrict *__restrict gauge , 
		   struct site *__restrict lat , 
		   const GLU_complex *__restrict *__restrict gauge_up ,
		   const size_t t )
{
  size_t i ;
  const size_t slice = LCU  *  t ; 
#pragma omp for private(i)
  PFOR(  i = 0  ;  i < LCU  ;  i ++ ) {
    const size_t j = slice + i ;
    #if ND == 4   
    inline_gtransform_local( gauge[i] , lat[j].O[0] , gauge[lat[i].neighbor[0]] ) ;
    inline_gtransform_local( gauge[i] , lat[j].O[1] , gauge[lat[i].neighbor[1]] ) ;
    inline_gtransform_local( gauge[i] , lat[j].O[2] , gauge[lat[i].neighbor[2]] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND - 1  ; mu++ ) {
      inline_gtransform_local( gauge[i] , lat[j].O[mu] , gauge[lat[i].neighbor[mu]] ) ;
    }
    #endif
    inline_gtransform_local( gauge[i] , lat[j].O[ND-1] , gauge_up[i] ) ;
  }
  return ;
}
