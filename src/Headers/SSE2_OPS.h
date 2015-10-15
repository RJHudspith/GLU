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
#define SSE2_MULCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_movedup_pd( a ) , b ) , \
					_mm_mul_pd( _mm_unpackhi_pd( a , SSE_FLIP(a) ) , \
						    _mm_shuffle_pd( b , b , 1 ) ) ) )
#else
#define SSE2_MULCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , a ) , b ) , \
					_mm_mul_pd( _mm_unpackhi_pd( a , SSE_FLIP(a) ) , \
						    _mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// performs a * b using intrinsics
#ifdef __SSE3__
#define SSE2_MUL(a,b) ( _mm_addsub_pd( _mm_mul_pd( _mm_movedup_pd( a ) , b ) , \
				       _mm_mul_pd( _mm_unpackhi_pd( a , a ) , \
						   _mm_shuffle_pd( b , b , 1 ) ) ) )
#else
#define SSE2_MUL(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , a ) , b ) , \
				    _mm_mul_pd( _mm_unpackhi_pd( SSE_FLIP(a) , a ) , \
						_mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// performs a * conj( b ) using intrinsics
#ifdef __SSE3__
#define SSE2_MUL_CONJ(a,b) ( _mm_addsub_pd( _mm_mul_pd( _mm_unpackhi_pd( a , a ) , \
							_mm_shuffle_pd( b , b , 1 ) ) , \
					    _mm_mul_pd( _mm_movedup_pd( a ) , SSE_FLIP( b ) ) ) )
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
#define SSE2_CONJ(a) ( _mm_move_sd( SSE_FLIP(a) , a ) ) ;

#endif
