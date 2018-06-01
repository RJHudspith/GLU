/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (SSE2_OPS.h) is part of GLU.

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
   @file SSE2_OPS.h
   @brief macro definitions for the SSE2 -> SSE4.2 operations we use
 */
#ifndef SSE2_OPS_H
#define SSE2_OPS_H

// gcc/clang allow for + / - * with SSE types, icc does not
#if !(defined __ICC)
  #define SSE_FLIP(a) ( -a )
#else
  #define SSE_FLIP(a) ( _mm_xor_pd( a , _mm_set1_pd( -0.0 ) ) )
#endif

// performs conj(a) * b using intrinsics
#ifdef __SSE3__
  #ifdef __FMA__
  #define SSE2_MULCONJ(a,b) ( _mm_fmadd_pd(  _mm_movedup_pd( a ) , b , \
					     _mm_mul_pd( _mm_unpackhi_pd( a , SSE_FLIP(a) ) , \
							 _mm_shuffle_pd( b , b , 1 ) ) ) )
  #else
  #define SSE2_MULCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_movedup_pd( a ) , b ) , \
					  _mm_mul_pd( _mm_unpackhi_pd( a , SSE_FLIP(a) ) , \
						      _mm_shuffle_pd( b , b , 1 ) ) ) )
  #endif
#else
#define SSE2_MULCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , a ) , b ) , \
					_mm_mul_pd( _mm_unpackhi_pd( a , SSE_FLIP(a) ) , \
						    _mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// performs a * b using intrinsics
#ifdef __SSE3__
  #ifdef __FMA__
  #define SSE2_MUL(a,b) ( _mm_fmaddsub_pd( _mm_movedup_pd( a ) , b , \
					   _mm_mul_pd( _mm_unpackhi_pd( a , a ) , \
						       _mm_shuffle_pd( b , b , 1 ) ) ) )
  #else
  #define SSE2_MUL(a,b) ( _mm_addsub_pd( _mm_mul_pd( _mm_movedup_pd( a ) , b ) , \
					 _mm_mul_pd( _mm_unpackhi_pd( a , a ) , \
						     _mm_shuffle_pd( b , b , 1 ) ) ) )
  #endif
#else
#define SSE2_MUL(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , a ) , b ) , \
				    _mm_mul_pd( _mm_unpackhi_pd( SSE_FLIP(a) , a ) , \
						_mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// performs a * conj( b ) using intrinsics
#ifdef __SSE3__
  #ifdef __FMA__
  #define SSE2_MUL_CONJ(a,b) ( _mm_fmaddsub_pd(  _mm_unpackhi_pd( a , a ) , \
						 _mm_shuffle_pd( b , b , 1 ) , \
						 _mm_mul_pd( _mm_movedup_pd( a ) , SSE_FLIP( b ) ) ) )
  #else
  #define SSE2_MUL_CONJ(a,b) ( _mm_addsub_pd( _mm_mul_pd( _mm_unpackhi_pd( a , a ) , \
							  _mm_shuffle_pd( b , b , 1 ) ) , \
					      _mm_mul_pd( _mm_movedup_pd( a ) , SSE_FLIP( b ) ) ) )
  #endif
#else
#define SSE2_MUL_CONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , SSE_FLIP(a) ) , b ) , \
					 _mm_mul_pd( _mm_unpackhi_pd( a , a ) , \
						     _mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// performs conj( a ) * conj( b ) using intrinsics
#ifdef __SSE3__
#define SSE2_MUL_CONJCONJ(a,b) ( _mm_hsub_pd( _mm_mul_pd( a , b ) , \
					      _mm_mul_pd( a , _mm_shuffle_pd( SSE_FLIP(b) , b , 1 ) ) ) )
#else
#define SSE2_MUL_CONJCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , SSE_FLIP(a) ) , b ) , \
					     _mm_mul_pd( _mm_unpackhi_pd( SSE_FLIP(a) , SSE_FLIP(a) ) , \
							 _mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// multiply by I
#define SSE2_iMUL(a) ( _mm_shuffle_pd( SSE_FLIP(a) , a , 1 ) )

// multiply by -I
#define SSE2_miMUL(a) ( _mm_shuffle_pd( a , SSE_FLIP(a) , 1 ) )

// complex conjugate
#define SSE2_CONJ(a) ( _mm_move_sd( SSE_FLIP(a) , a ) )

#endif
