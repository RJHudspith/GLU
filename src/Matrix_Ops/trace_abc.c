/**
   @file traces_abc.c
   @brief trace of the product of three matrices
 */
#include "Mainfile.h"

#if !(defined HAVE_IMMINTRIN_H) || (defined SINGLE_PREC)

double
Re_trace_abc_dag_suNC( const GLU_complex a[ NCNC ] , 
		       const GLU_complex b[ NCNC ] , 
		       const GLU_complex c[ NCNC ] )
{
  register double sum ;
#if NC == 3
  //const GLU_complex a0 = ( a[0] * b[0] + a[1] * b[3] + a[2] * b[6] ) ;
  const GLU_real Ra0 = creal( a[0] ) * creal( b[0] ) - cimag( a[0] ) * cimag( b[0] ) +\
                       creal( a[1] ) * creal( b[3] ) - cimag( a[1] ) * cimag( b[3] ) +\
                       creal( a[2] ) * creal( b[6] ) - cimag( a[2] ) * cimag( b[6] ) ;
  const GLU_real Ia0 = creal( a[0] ) * cimag( b[0] ) + cimag( a[0] ) * creal( b[0] ) +\
                       creal( a[1] ) * cimag( b[3] ) + cimag( a[1] ) * creal( b[3] ) +\
                       creal( a[2] ) * cimag( b[6] ) + cimag( a[2] ) * creal( b[6] ) ;
  //const GLU_complex a1 = ( a[0] * b[1] + a[1] * b[4] + a[2] * b[7] ) ;
  const GLU_real Ra1 = creal( a[0] ) * creal( b[1] ) - cimag( a[0] ) * cimag( b[1] ) +\
                       creal( a[1] ) * creal( b[4] ) - cimag( a[1] ) * cimag( b[4] ) +\
                       creal( a[2] ) * creal( b[7] ) - cimag( a[2] ) * cimag( b[7] ) ;
  const GLU_real Ia1 = creal( a[0] ) * cimag( b[1] ) + cimag( a[0] ) * creal( b[1] ) +\
                       creal( a[1] ) * cimag( b[4] ) + cimag( a[1] ) * creal( b[4] ) +\
                       creal( a[2] ) * cimag( b[7] ) + cimag( a[2] ) * creal( b[7] ) ;
  //const GLU_complex a2 = ( a[0] * b[2] + a[1] * b[5] + a[2] * b[8] ) ;
  const GLU_real Ra2 = creal( a[0] ) * creal( b[2] ) - cimag( a[0] ) * cimag( b[2] ) +\
                       creal( a[1] ) * creal( b[5] ) - cimag( a[1] ) * cimag( b[5] ) +\
                       creal( a[2] ) * creal( b[8] ) - cimag( a[2] ) * cimag( b[8] ) ;
  const GLU_real Ia2 = creal( a[0] ) * cimag( b[2] ) + cimag( a[0] ) * creal( b[2] ) +\
                       creal( a[1] ) * cimag( b[5] ) + cimag( a[1] ) * creal( b[5] ) +\
                       creal( a[2] ) * cimag( b[8] ) + cimag( a[2] ) * creal( b[8] ) ;
  //const GLU_complex a3 = ( a[3] * b[0] + a[4] * b[3] + a[5] * b[6] ) ;
  const GLU_real Ra3 = creal( a[3] ) * creal( b[0] ) - cimag( a[3] ) * cimag( b[0] ) +\
                       creal( a[4] ) * creal( b[3] ) - cimag( a[4] ) * cimag( b[3] ) +\
                       creal( a[5] ) * creal( b[6] ) - cimag( a[5] ) * cimag( b[6] ) ;
  const GLU_real Ia3 = creal( a[3] ) * cimag( b[0] ) + cimag( a[3] ) * creal( b[0] ) +\
                       creal( a[4] ) * cimag( b[3] ) + cimag( a[4] ) * creal( b[3] ) +\
                       creal( a[5] ) * cimag( b[6] ) + cimag( a[5] ) * creal( b[6] ) ;
  //const GLU_complex a4 = ( a[3] * b[1] + a[4] * b[4] + a[5] * b[7] ) ;
  const GLU_real Ra4 = creal( a[3] ) * creal( b[1] ) - cimag( a[3] ) * cimag( b[1] ) +\
                       creal( a[4] ) * creal( b[4] ) - cimag( a[4] ) * cimag( b[4] ) +\
                       creal( a[5] ) * creal( b[7] ) - cimag( a[5] ) * cimag( b[7] ) ;
  const GLU_real Ia4 = creal( a[3] ) * cimag( b[1] ) + cimag( a[3] ) * creal( b[1] ) +\
                       creal( a[4] ) * cimag( b[4] ) + cimag( a[4] ) * creal( b[4] ) +\
                       creal( a[5] ) * cimag( b[7] ) + cimag( a[5] ) * creal( b[7] ) ;
  //const GLU_complex a5 = ( a[3] * b[2] + a[4] * b[5] + a[5] * b[8] ) ;
  const GLU_real Ra5 = creal( a[3] ) * creal( b[2] ) - cimag( a[3] ) * cimag( b[2] ) +\
                       creal( a[4] ) * creal( b[5] ) - cimag( a[4] ) * cimag( b[5] ) +\
                       creal( a[5] ) * creal( b[8] ) - cimag( a[5] ) * cimag( b[8] ) ;
  const GLU_real Ia5 = creal( a[3] ) * cimag( b[2] ) + cimag( a[3] ) * creal( b[2] ) +\
                       creal( a[4] ) * cimag( b[5] ) + cimag( a[4] ) * creal( b[5] ) +\
                       creal( a[5] ) * cimag( b[8] ) + cimag( a[5] ) * creal( b[8] ) ;
  //const GLU_complex a6 = conj( a1 * a5 - a2 * a4 ) ;
  const GLU_real Ra6 = +( Ra1 * Ra5 - Ia1 * Ia5 - Ra2 * Ra4 + Ia2 * Ia4 ) ;
  const GLU_real Ia6 = -( Ra1 * Ia5 + Ia1 * Ra5 - Ra2 * Ia4 - Ia2 * Ra4 ) ;
  //const GLU_complex a7 = conj( a2 * a3 - a0 * a5 ) ;
  const GLU_real Ra7 = +( Ra2 * Ra3 - Ia2 * Ia3 - Ra0 * Ra5 + Ia0 * Ia5 ) ;
  const GLU_real Ia7 = -( Ra2 * Ia3 + Ia2 * Ra3 - Ra0 * Ia5 - Ia0 * Ra5 ) ;
  //const GLU_complex a8 = conj( a0 * a4 - a1 * a3 ) ;
  const GLU_real Ra8 = +( Ra0 * Ra4 - Ia0 * Ia4 - Ra1 * Ra3 + Ia1 * Ia3 ) ;
  const GLU_real Ia8 = -( Ra0 * Ia4 + Ia0 * Ra4 - Ra1 * Ia3 - Ia1 * Ra3 ) ;
  // and compute the trace
  sum = Ra0 * creal( c[0] ) + Ia0 * cimag( c[0] ) ; 
  sum = Ra1 * creal( c[1] ) + Ia1 * cimag( c[1] ) + sum ; 
  sum = Ra2 * creal( c[2] ) + Ia2 * cimag( c[2] ) + sum ; 
  sum = Ra3 * creal( c[3] ) + Ia3 * cimag( c[3] ) + sum ; 
  sum = Ra4 * creal( c[4] ) + Ia4 * cimag( c[4] ) + sum ; 
  sum = Ra5 * creal( c[5] ) + Ia5 * cimag( c[5] ) + sum ; 
  sum = Ra6 * creal( c[6] ) + Ia6 * cimag( c[6] ) + sum ; 
  sum = Ra7 * creal( c[7] ) + Ia7 * cimag( c[7] ) + sum ; 
  sum = Ra8 * creal( c[8] ) + Ia8 * cimag( c[8] ) + sum ; 
#elif NC == 2
  const GLU_complex a0 = a[0] * b[0] + a[1] * b[2] ;
  const GLU_complex a1 = a[0] * b[1] + a[1] * b[3] ;
  sum  = creal( a0 ) * creal( c[0] ) + cimag( a0 ) * cimag( c[0] ) ;
  sum += creal( a1 ) * creal( c[1] ) + cimag( a1 ) * cimag( c[1] ) ;
  sum -= creal( a1 ) * creal( c[2] ) - cimag( a1 ) * cimag( c[2] ) ;
  sum += creal( a0 ) * creal( c[3] ) - cimag( a0 ) * cimag( c[3] ) ;
#else
  GLU_real trABCdag ;
  trace_abc_dag_Re( &trABCdag , a , b , c ) ;
  sum = (double)trABCdag ;
#endif
  return sum ;
}

// Trace of the product of three matrices //
void
trace_abc( GLU_complex *__restrict tr , 
	   const GLU_complex a[ NCNC ] , 
	   const GLU_complex b[ NCNC ] , 
	   const GLU_complex c[ NCNC ] )
{
#if NC == 3
  *tr = ( a[0] * b[0] + a[1] * b[3] + a[2] * b[6] ) * c[0] +\
    ( a[3] * b[0] + a[4] * b[3] + a[5] * b[6] ) * c[1] +    \
    ( a[6] * b[0] + a[7] * b[3] + a[8] * b[6] ) * c[2] +    \
    ( a[0] * b[1] + a[1] * b[4] + a[2] * b[7] ) * c[3] +    \
    ( a[3] * b[1] + a[4] * b[4] + a[5] * b[7] ) * c[4] +    \
    ( a[6] * b[1] + a[7] * b[4] + a[8] * b[7] ) * c[5] +    \
    ( a[0] * b[2] + a[1] * b[5] + a[2] * b[8] ) * c[6] +    \
    ( a[3] * b[2] + a[4] * b[5] + a[5] * b[8] ) * c[7] +    \
    ( a[6] * b[2] + a[7] * b[5] + a[8] * b[8] ) * c[8] ;
#elif NC == 2
  *tr = ( a[0] * b[0] + a[1] * b[2] ) * c[0] +	\
    ( a[2] * b[0] + a[3] * b[2] ) * c[1] +	\
    ( a[0] * b[1] + a[1] * b[3] ) * c[2] +	\
    ( a[2] * b[1] + a[3] * b[3] ) * c[3] ;
#else
  const GLU_complex *pB , *pA ;
  register GLU_real sumr = 0.0 , sumi = 0.0 ;
  GLU_real insumr = 0.0 , insumi = 0.0 ;
  size_t i , j , k ;
  for( i = 0 ; i < NC ; i++ ) {
    pA = a ;
    for( j = 0 ; j < NC ; j++ ) {
      pB = b ;
      insumr = insumi = 0.0 ;
      for( k = 0 ; k < NC ; k++ ) {
	// unroll the mul
	insumr += ( creal( pA[k] ) * creal( pB[i] ) -
		    cimag( pA[k] ) * cimag( pB[i] ) ) ; 
	insumi += ( creal( pA[k] ) * cimag( pB[i] ) + 
		    cimag( pA[k] ) * creal( pB[i] ) ) ;
	pB += NC ;
      }
      sumr += insumr * creal( c[j+i*NC] ) - insumi * cimag( c[j+i*NC] ) ;
      sumi += insumr * cimag( c[j+i*NC] ) + insumi * creal( c[j+i*NC] ) ;
      pA += NC ;
    }
  }
  *tr = sumr + I * sumi ;
#endif
  return ;
}

// this is trace( a . b . c^{\dagger} )
void
trace_abc_dag( GLU_complex *__restrict tr , 
	       const GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC ] , 
	       const GLU_complex c[ NCNC ] )
{
#if NC == 3
  *tr = ( a[0] * b[0] + a[1] * b[3] + a[2] * b[6] ) * conj( c[0] ) +	\
    ( a[3] * b[0] + a[4] * b[3] + a[5] * b[6] ) * conj( c[3] ) +	\
    ( a[6] * b[0] + a[7] * b[3] + a[8] * b[6] ) * conj( c[6] ) +	\
    ( a[0] * b[1] + a[1] * b[4] + a[2] * b[7] ) * conj( c[1] ) +	\
    ( a[3] * b[1] + a[4] * b[4] + a[5] * b[7] ) * conj( c[4] ) +	\
    ( a[6] * b[1] + a[7] * b[4] + a[8] * b[7] ) * conj( c[7] ) +	\
    ( a[0] * b[2] + a[1] * b[5] + a[2] * b[8] ) * conj( c[2] ) +	\
    ( a[3] * b[2] + a[4] * b[5] + a[5] * b[8] ) * conj( c[5] ) +	\
    ( a[6] * b[2] + a[7] * b[5] + a[8] * b[8] ) * conj( c[8] ) ;
#elif NC == 2
  *tr = ( a[0] * b[0] + a[1] * b[2] ) * conj( c[0] ) + \
    ( a[0] * b[1] + a[1] * b[3] ) * conj( c[1] ) + \
    ( a[2] * b[0] + a[3] * b[2] ) * conj( c[2] ) + \
    ( a[2] * b[1] + a[3] * b[3] ) * conj( c[3] ) ;
#else
  const GLU_complex *pB , *pA ;
  register GLU_real sumr = 0.0 , sumi = 0.0 ;
  GLU_real insumr = 0.0 , insumi = 0.0 ;
  size_t i , j , k ;
  for( i = 0 ; i < NC ; i++ ) {
    pA = a ;
    for( j = 0 ; j < NC ; j++ ) {
      pB = b ;
      insumr = insumi = 0.0 ;
      for( k = 0 ; k < NC ; k++ ) {
	// unroll the mul
	insumr += ( creal( pA[k] ) * creal( pB[i] ) -
		    cimag( pA[k] ) * cimag( pB[i] ) ) ; 
	insumi += ( creal( pA[k] ) * cimag( pB[i] ) + 
		    cimag( pA[k] ) * creal( pB[i] ) ) ;
	pB += NC ;
      }
      sumr +=  insumr * creal( c[i+j*NC] ) + insumi * cimag( c[i+j*NC] ) ;
      sumi += -insumr * cimag( c[i+j*NC] ) + insumi * creal( c[i+j*NC] ) ;
      pA += NC ;
    }
  }
  *tr = sumr + I * sumi ;
#endif
  return ;
}

// this is trace( a . b . c^{\dagger} )
void
trace_abc_dag_Re( GLU_real *__restrict tr , 
		  const GLU_complex a[ NCNC ] , 
		  const GLU_complex b[ NCNC ] , 
		  const GLU_complex c[ NCNC ] )
{
#if NC == 3
  *tr = creal( ( a[0] * b[0] + a[1] * b[3] + a[2] * b[6] ) * conj( c[0] ) + \
	       ( a[0] * b[1] + a[1] * b[4] + a[2] * b[7] ) * conj( c[1] ) + \
	       ( a[0] * b[2] + a[1] * b[5] + a[2] * b[8] ) * conj( c[2] ) + \
	       ( a[3] * b[0] + a[4] * b[3] + a[5] * b[6] ) * conj( c[3] ) + \
	       ( a[3] * b[1] + a[4] * b[4] + a[5] * b[7] ) * conj( c[4] ) + \
	       ( a[3] * b[2] + a[4] * b[5] + a[5] * b[8] ) * conj( c[5] ) + \
	       ( a[6] * b[0] + a[7] * b[3] + a[8] * b[6] ) * conj( c[6] ) + \
	       ( a[6] * b[1] + a[7] * b[4] + a[8] * b[7] ) * conj( c[7] ) + \
	       ( a[6] * b[2] + a[7] * b[5] + a[8] * b[8] ) * conj( c[8] ) ) ;
#elif NC == 2
  *tr = creal( ( a[0] * b[0] + a[1] * b[2] ) * conj( c[0] ) +	\
	       ( a[0] * b[1] + a[1] * b[3] ) * conj( c[1] ) +	\
	       ( a[2] * b[0] + a[3] * b[2] ) * conj( c[2] ) +	\
	       ( a[2] * b[1] + a[3] * b[3] ) * conj( c[3] ) ) ;
#else
  const GLU_complex *pB , *pA ;
  register double sumr = 0.0 ;
  double insumr = 0.0 , insumi = 0.0 ;
  size_t i , j , k ;
  for( i = 0 ; i < NC ; i++ ) {
    pA = a ;
    for( j = 0 ; j < NC ; j++ ) {
      pB = b ;
      insumr = insumi = 0.0 ;
      for( k = 0 ; k < NC ; k++ ) {
	// unroll the mul
	insumr += ( creal( pA[k] ) * creal( pB[i] ) -
		    cimag( pA[k] ) * cimag( pB[i] ) ) ; 
	insumi += ( creal( pA[k] ) * cimag( pB[i] ) + 
		    cimag( pA[k] ) * creal( pB[i] ) ) ;
	pB += NC ;
      }
      sumr +=  insumr * creal( c[i+j*NC] ) + insumi * cimag( c[i+j*NC] ) ;
      pA += NC ;
    }
  }
  *tr = sumr ;
#endif
  return ;
}

#endif
