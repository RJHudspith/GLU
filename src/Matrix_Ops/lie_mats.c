/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (lie_mats.c) is part of GLU.

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
   @file lie_mats.c
   @brief has the lie-matrix definitions and structure constants for fundamental SU(N)
 */

#include "Mainfile.h"

// general structure function
struct struc_func {
  size_t a, b , c ;
  GLU_real val ;
} ;

// f^{abc} and d^{abc}, get allocated or not
static struct struc_func *f = NULL ; 
static struct struc_func *d = NULL ; 

// the generators get allocated at some point or not
static GLU_complex **lambda = NULL ; 

// some counters
static size_t fcount = 0 , dcount = 0 ;

// inline for the matrix idx
static inline int mat_idx( const size_t row , const size_t column ) { return row + column * NC ; }

#if NC > 3
// Obviously,
// Tr[ ( A.B - B.A ).C ] == Tr[ A.B.C ] - Tr[ B.A.C ]
// and then we can use our solution for the trace of three arbitrary
// matrices to remove a bunch of ops.
static void
compute_f_and_d( GLU_real *f , 
		 GLU_real *d ,
		 const GLU_complex lambdaA[ NCNC ] ,
		 const GLU_complex lambdaB[ NCNC ] ,
		 const GLU_complex lambdaC[ NCNC ] )
{
  GLU_complex tr1 , tr2 ;
  trace_abc( &tr1 , lambdaA , lambdaB , lambdaC ) ;
  trace_abc( &tr2 , lambdaB , lambdaA , lambdaC ) ;
  *f = -2.0 * I * ( tr1 - tr2 ) ;
  *d = 2.0 * ( tr1 + tr2 ) ;
  return ;
}
#endif

static void
actually_compute_fs_and_ds( void )
{
  // I did not hand type this, still. Yuck!
#if NC == 3
  fcount = 54 ;
  dcount = 58 ;
  f = ( struct struc_func* )malloc( fcount * sizeof( struct struc_func ) ) ;
  d = ( struct struc_func* )malloc( dcount * sizeof( struct struc_func ) ) ;
  const register GLU_real r3_two = 0.8660254037844386 ;
  f[0].val = 1.0      ; f[0].a = 0 ; f[0].b = 1 ; f[0].c = 2 ; 
  f[1].val = -1.0     ; f[1].a = 0 ; f[1].b = 2 ; f[1].c = 1 ; 
  f[2].val = 0.5      ; f[2].a = 0 ; f[2].b = 3 ; f[2].c = 6 ; 
  f[3].val = -0.5     ; f[3].a = 0 ; f[3].b = 4 ; f[3].c = 5 ; 
  f[4].val = 0.5      ; f[4].a = 0 ; f[4].b = 5 ; f[4].c = 4 ; 
  f[5].val = -0.5     ; f[5].a = 0 ; f[5].b = 6 ; f[5].c = 3 ; 
  f[6].val = -1.0     ; f[6].a = 1 ; f[6].b = 0 ; f[6].c = 2 ; 
  f[7].val = 1.0      ; f[7].a = 1 ; f[7].b = 2 ; f[7].c = 0 ; 
  f[8].val = 0.5      ; f[8].a = 1 ; f[8].b = 3 ; f[8].c = 5 ; 
  f[9].val = 0.5      ; f[9].a = 1 ; f[9].b = 4 ; f[9].c = 6 ; 
  f[10].val = -0.5    ; f[10].a = 1 ; f[10].b = 5 ; f[10].c = 3 ; 
  f[11].val = -0.5    ; f[11].a = 1 ; f[11].b = 6 ; f[11].c = 4 ; 
  f[12].val = 1.0     ; f[12].a = 2 ; f[12].b = 0 ; f[12].c = 1 ; 
  f[13].val = -1.0    ; f[13].a = 2 ; f[13].b = 1 ; f[13].c = 0 ; 
  f[14].val = 0.5     ; f[14].a = 2 ; f[14].b = 3 ; f[14].c = 4 ; 
  f[15].val = -0.5    ; f[15].a = 2 ; f[15].b = 4 ; f[15].c = 3 ; 
  f[16].val = -0.5    ; f[16].a = 2 ; f[16].b = 5 ; f[16].c = 6 ; 
  f[17].val = 0.5     ; f[17].a = 2 ; f[17].b = 6 ; f[17].c = 5 ; 
  f[18].val = -0.5    ; f[18].a = 3 ; f[18].b = 0 ; f[18].c = 6 ; 
  f[19].val = -0.5    ; f[19].a = 3 ; f[19].b = 1 ; f[19].c = 5 ; 
  f[20].val = -0.5    ; f[20].a = 3 ; f[20].b = 2 ; f[20].c = 4 ; 
  f[21].val = 0.5     ; f[21].a = 3 ; f[21].b = 4 ; f[21].c = 2 ; 
  f[22].val = r3_two  ; f[22].a = 3 ; f[22].b = 4 ; f[22].c = 7 ; 
  f[23].val = 0.5     ; f[23].a = 3 ; f[23].b = 5 ; f[23].c = 1 ; 
  f[24].val = 0.5     ; f[24].a = 3 ; f[24].b = 6 ; f[24].c = 0 ; 
  f[25].val = -r3_two ; f[25].a = 3 ; f[25].b = 7 ; f[25].c = 4 ; 
  f[26].val = 0.5     ; f[26].a = 4 ; f[26].b = 0 ; f[26].c = 5 ; 
  f[27].val = -0.5    ; f[27].a = 4 ; f[27].b = 1 ; f[27].c = 6 ; 
  f[28].val = 0.5     ; f[28].a = 4 ; f[28].b = 2 ; f[28].c = 3 ; 
  f[29].val = -0.5    ; f[29].a = 4 ; f[29].b = 3 ; f[29].c = 2 ; 
  f[30].val = -r3_two ; f[30].a = 4 ; f[30].b = 3 ; f[30].c = 7 ; 
  f[31].val = -0.5    ; f[31].a = 4 ; f[31].b = 5 ; f[31].c = 0 ; 
  f[32].val = 0.5     ; f[32].a = 4 ; f[32].b = 6 ; f[32].c = 1 ; 
  f[33].val = r3_two  ; f[33].a = 4 ; f[33].b = 7 ; f[33].c = 3 ; 
  f[34].val = -0.5    ; f[34].a = 5 ; f[34].b = 0 ; f[34].c = 4 ; 
  f[35].val = 0.5     ; f[35].a = 5 ; f[35].b = 1 ; f[35].c = 3 ; 
  f[36].val = 0.5     ; f[36].a = 5 ; f[36].b = 2 ; f[36].c = 6 ; 
  f[37].val = -0.5    ; f[37].a = 5 ; f[37].b = 3 ; f[37].c = 1 ; 
  f[38].val = 0.5     ; f[38].a = 5 ; f[38].b = 4 ; f[38].c = 0 ; 
  f[39].val = -0.5    ; f[39].a = 5 ; f[39].b = 6 ; f[39].c = 2 ; 
  f[40].val = r3_two  ; f[40].a = 5 ; f[40].b = 6 ; f[40].c = 7 ; 
  f[41].val = -r3_two ; f[41].a = 5 ; f[41].b = 7 ; f[41].c = 6 ; 
  f[42].val = 0.5     ; f[42].a = 6 ; f[42].b = 0 ; f[42].c = 3 ; 
  f[43].val = 0.5     ; f[43].a = 6 ; f[43].b = 1 ; f[43].c = 4 ; 
  f[44].val = -0.5    ; f[44].a = 6 ; f[44].b = 2 ; f[44].c = 5 ; 
  f[45].val = -0.5    ; f[45].a = 6 ; f[45].b = 3 ; f[45].c = 0 ; 
  f[46].val = -0.5    ; f[46].a = 6 ; f[46].b = 4 ; f[46].c = 1 ; 
  f[47].val = 0.5     ; f[47].a = 6 ; f[47].b = 5 ; f[47].c = 2 ; 
  f[48].val = -r3_two ; f[48].a = 6 ; f[48].b = 5 ; f[48].c = 7 ; 
  f[49].val = r3_two  ; f[49].a = 6 ; f[49].b = 7 ; f[49].c = 5 ; 
  f[50].val = r3_two  ; f[50].a = 7 ; f[50].b = 3 ; f[50].c = 4 ; 
  f[51].val = -r3_two ; f[51].a = 7 ; f[51].b = 4 ; f[51].c = 3 ; 
  f[52].val = r3_two  ; f[52].a = 7 ; f[52].b = 5 ; f[52].c = 6 ; 
  f[53].val = -r3_two ; f[53].a = 7 ; f[53].b = 6 ; f[53].c = 5 ; 
  // and the d's
  register const GLU_real one_r3 = 0.57735026918962584 ;
  register const GLU_real one_twor3 = 0.28867513459481292 ;
  d[0].val = one_r3      ; d[0].a = 0 ; d[0].b = 0 ; d[0].c = 7 ; 
  d[1].val = 0.5         ; d[1].a = 0 ; d[1].b = 3 ; d[1].c = 5 ; 
  d[2].val = 0.5         ; d[2].a = 0 ; d[2].b = 4 ; d[2].c = 6 ; 
  d[3].val = 0.5         ; d[3].a = 0 ; d[3].b = 5 ; d[3].c = 3 ; 
  d[4].val = 0.5         ; d[4].a = 0 ; d[4].b = 6 ; d[4].c = 4 ; 
  d[5].val = one_r3      ; d[5].a = 0 ; d[5].b = 7 ; d[5].c = 0 ; 
  d[6].val = one_r3      ; d[6].a = 1 ; d[6].b = 1 ; d[6].c = 7 ; 
  d[7].val = -0.5        ; d[7].a = 1 ; d[7].b = 3 ; d[7].c = 6 ; 
  d[8].val = 0.5         ; d[8].a = 1 ; d[8].b = 4 ; d[8].c = 5 ; 
  d[9].val = 0.5         ; d[9].a = 1 ; d[9].b = 5 ; d[9].c = 4 ; 
  d[10].val = -0.5       ; d[10].a = 1 ; d[10].b = 6 ; d[10].c = 3 ; 
  d[11].val = one_r3     ; d[11].a = 1 ; d[11].b = 7 ; d[11].c = 1 ; 
  d[12].val = one_r3     ; d[12].a = 2 ; d[12].b = 2 ; d[12].c = 7 ; 
  d[13].val = 0.5        ; d[13].a = 2 ; d[13].b = 3 ; d[13].c = 3 ; 
  d[14].val = 0.5        ; d[14].a = 2 ; d[14].b = 4 ; d[14].c = 4 ; 
  d[15].val = -0.5       ; d[15].a = 2 ; d[15].b = 5 ; d[15].c = 5 ; 
  d[16].val = -0.5       ; d[16].a = 2 ; d[16].b = 6 ; d[16].c = 6 ; 
  d[17].val = one_r3     ; d[17].a = 2 ; d[17].b = 7 ; d[17].c = 2 ; 
  d[18].val = 0.5        ; d[18].a = 3 ; d[18].b = 0 ; d[18].c = 5 ; 
  d[19].val = -0.5       ; d[19].a = 3 ; d[19].b = 1 ; d[19].c = 6 ; 
  d[20].val = 0.5        ; d[20].a = 3 ; d[20].b = 2 ; d[20].c = 3 ; 
  d[21].val = 0.5        ; d[21].a = 3 ; d[21].b = 3 ; d[21].c = 2 ; 
  d[22].val = -one_twor3 ; d[22].a = 3 ; d[22].b = 3 ; d[22].c = 7 ; 
  d[23].val = 0.5        ; d[23].a = 3 ; d[23].b = 5 ; d[23].c = 0 ; 
  d[24].val = -0.5       ; d[24].a = 3 ; d[24].b = 6 ; d[24].c = 1 ; 
  d[25].val = -one_twor3 ; d[25].a = 3 ; d[25].b = 7 ; d[25].c = 3 ; 
  d[26].val = 0.5        ; d[26].a = 4 ; d[26].b = 0 ; d[26].c = 6 ; 
  d[27].val = 0.5        ; d[27].a = 4 ; d[27].b = 1 ; d[27].c = 5 ; 
  d[28].val = 0.5        ; d[28].a = 4 ; d[28].b = 2 ; d[28].c = 4 ; 
  d[29].val = 0.5        ; d[29].a = 4 ; d[29].b = 4 ; d[29].c = 2 ; 
  d[30].val = -one_twor3 ; d[30].a = 4 ; d[30].b = 4 ; d[30].c = 7 ; 
  d[31].val = 0.5        ; d[31].a = 4 ; d[31].b = 5 ; d[31].c = 1 ; 
  d[32].val = 0.5        ; d[32].a = 4 ; d[32].b = 6 ; d[32].c = 0 ; 
  d[33].val = -one_twor3 ; d[33].a = 4 ; d[33].b = 7 ; d[33].c = 4 ; 
  d[34].val = 0.5        ; d[34].a = 5 ; d[34].b = 0 ; d[34].c = 3 ; 
  d[35].val = 0.5        ; d[35].a = 5 ; d[35].b = 1 ; d[35].c = 4 ; 
  d[36].val = -0.5       ; d[36].a = 5 ; d[36].b = 2 ; d[36].c = 5 ; 
  d[37].val = 0.5        ; d[37].a = 5 ; d[37].b = 3 ; d[37].c = 0 ; 
  d[38].val = 0.5        ; d[38].a = 5 ; d[38].b = 4 ; d[38].c = 1 ; 
  d[39].val = -0.5       ; d[39].a = 5 ; d[39].b = 5 ; d[39].c = 2 ; 
  d[40].val = -one_twor3 ; d[40].a = 5 ; d[40].b = 5 ; d[40].c = 7 ; 
  d[41].val = -one_twor3 ; d[41].a = 5 ; d[41].b = 7 ; d[41].c = 5 ; 
  d[42].val = 0.5        ; d[42].a = 6 ; d[42].b = 0 ; d[42].c = 4 ; 
  d[43].val = -0.5       ; d[43].a = 6 ; d[43].b = 1 ; d[43].c = 3 ; 
  d[44].val = -0.5       ; d[44].a = 6 ; d[44].b = 2 ; d[44].c = 6 ; 
  d[45].val = -0.5       ; d[45].a = 6 ; d[45].b = 3 ; d[45].c = 1 ; 
  d[46].val = 0.5        ; d[46].a = 6 ; d[46].b = 4 ; d[46].c = 0 ; 
  d[47].val = -0.5       ; d[47].a = 6 ; d[47].b = 6 ; d[47].c = 2 ; 
  d[48].val = -one_twor3 ; d[48].a = 6 ; d[48].b = 6 ; d[48].c = 7 ; 
  d[49].val = -one_twor3 ; d[49].a = 6 ; d[49].b = 7 ; d[49].c = 6 ; 
  d[50].val = one_r3     ; d[50].a = 7 ; d[50].b = 0 ; d[50].c = 0 ; 
  d[51].val = one_r3     ; d[51].a = 7 ; d[51].b = 1 ; d[51].c = 1 ; 
  d[52].val = one_r3     ; d[52].a = 7 ; d[52].b = 2 ; d[52].c = 2 ; 
  d[53].val = -one_twor3 ; d[53].a = 7 ; d[53].b = 3 ; d[53].c = 3 ; 
  d[54].val = -one_twor3 ; d[54].a = 7 ; d[54].b = 4 ; d[54].c = 4 ; 
  d[55].val = -one_twor3 ; d[55].a = 7 ; d[55].b = 5 ; d[55].c = 5 ; 
  d[56].val = -one_twor3 ; d[56].a = 7 ; d[56].b = 6 ; d[56].c = 6 ; 
  d[57].val = -one_r3    ; d[57].a = 7 ; d[57].b = 7 ; d[57].c = 7 ;
  // well that was gratuitous
#elif NC == 2
  fcount = 6 ;
  f = ( struct struc_func* )malloc( fcount * sizeof( struct struc_func ) ) ;
  f[0].val = 1.0  ; f[0].a = 0 ; f[0].b = 1 ; f[0].c = 2 ; 
  f[1].val = -1.0 ; f[1].a = 0 ; f[1].b = 2 ; f[1].c = 1 ; 
  f[2].val = -1.0 ; f[2].a = 1 ; f[2].b = 0 ; f[2].c = 2 ; 
  f[3].val = 1.0  ; f[3].a = 1 ; f[3].b = 2 ; f[3].c = 0 ; 
  f[4].val = 1.0  ; f[4].a = 2 ; f[4].b = 0 ; f[4].c = 1 ; 
  f[5].val = -1.0 ; f[5].a = 2 ; f[5].b = 1 ; f[5].c = 0 ;
#else
  size_t a , b , c ; 
  //  Two pass variant ... The first pass is to calculate how many there are and then the second to allocate
  //  this was a lot uglier than my previous idea of using a linked list, but this can be used in parallel
  //  easier than traversing a linked list when calling for the f's and d's, not saying it isn't possible
  GLU_real ff , dd ;
  for( a = 0 ; a < NCNC - 1 ; a++ ) {
    for( b = 0 ; b < NCNC - 1 ; b++ ) {
      for( c = 0 ; c < NCNC - 1 ; c++ ) {
	compute_f_and_d( &ff , &dd , lambda[a] , lambda[b] , lambda[c] ) ;
	if( fabs( dd ) > PREC_TOL ) { dcount++ ; }
	if( fabs( ff ) > PREC_TOL ) { fcount++ ; }
      }
    }
  }
  f = ( struct struc_func* )malloc( fcount * sizeof( struct struc_func ) ) ;
  d = ( struct struc_func* )malloc( dcount * sizeof( struct struc_func ) ) ;
  fcount = dcount = 0 ;
  for( a = 0 ; a < NCNC - 1 ; a++ ) {
    for( b = 0 ; b < NCNC - 1 ; b++ ) {
      for( c = 0 ; c < NCNC - 1 ; c++ ) {
	// compute both the f and the d
	compute_f_and_d( &ff , &dd , lambda[a] , lambda[b] , lambda[c] ) ;
	if( fabs( dd ) > PREC_TOL ) { 
	  d[ dcount ].val = dd ;
	  d[ dcount ].a = a ; d[ dcount ].b = b ; d[ dcount ].c = c ;
	  dcount ++ ;
	}
	if( fabs( ff ) > PREC_TOL ) { 
	  f[ fcount ].val = ff ;
	  f[ fcount ].a = a ; f[ fcount ].b = b ; f[ fcount ].c = c ;
	  fcount++ ; 
	}
	// end of the cloop
      }
    }
  }
#endif
  // rad, we are done.
  return ;
}

// computes the lie-elements
static INLINE_VOID
lie_data( GLU_complex a[ NCNC - 1 ]  ,
	  const GLU_complex A[ NCNC ] )
{
#if NC == 3
  a[0] = A[1] + A[3] ;
  a[1] = I * ( A[1] - A[3] ) ;
  a[2] = A[0] - A[4] ;
  a[3] = A[2] + A[6] ;
  a[4] = I * ( A[2] - A[6] ) ;
  a[5] = A[5] + A[7] ;
  a[6] = I * ( A[5] - A[7] ) ;
  a[7] = -A[8] * 1.7320508075688772 ; // multiply by sqrt 3
#elif NC == 2
  a[0] = 0.5 * ( A[1] + A[2] ) ;
  a[1] = 0.5 * ( I * ( A[1] - A[2] ) ) ;
  a[2] = ( A[0] ) ;
#else
  // calculation of this is not so bad, a[i] = 2 Tr( A T^{i} )
  size_t i ;
  for( i = 0 ; i < NCNC-1 ; i++ ) {
    GLU_complex tr ;
    trace_ab( &tr , A , lambda[i] ) ; 
    a[i] = 2.0 * (GLU_complex)tr ;
  }
#endif
  return ;
}

// compute them by looking them up or by actually computing them
void
compute_fs_and_ds( void )
{
  // precomputation look up
#ifdef NOT_CONDOR_MODE
  // have a look for a file
  char str[ 256 ] ;
  sprintf( str , "%s/Local/Moments/SU%d_fabc_dabc.config" , HAVE_PREFIX , NC ) ;
  FILE *fabc_dabc = fopen( str , "rb" ) ;
  fprintf( stdout , "[LIE] Checking for precomputed file %s ... \n" , str ) ;
  // first element is the length of f, second the length of d
  if( fabc_dabc == NULL ) {
    fprintf( stdout , "[LIE] File not found. Writing one.\n" ) ; 
    actually_compute_fs_and_ds( ) ;
    // write them out
    fabc_dabc = fopen( str , "wb" ) ;
    int temp[1] = { fcount } ;
    fwrite( temp , sizeof( int ) , 1 , fabc_dabc ) ;
    temp[0] = dcount ;
    fwrite( temp , sizeof( int ) , 1 , fabc_dabc ) ;
    fwrite( f , sizeof( struct struc_func ) , fcount , fabc_dabc ) ;
    fwrite( d , sizeof( struct struc_func ) , dcount , fabc_dabc ) ;
  } else {
    fprintf( stdout , "[LIE] File found (%s) ... Reading it \n" , str ) ; 
    if( fread( &fcount , sizeof( int ) , 1 , fabc_dabc ) != 1 ) return ;
    if( fread( &dcount , sizeof( int ) , 1 , fabc_dabc ) != 1 ) return ;
    f = ( struct struc_func* )malloc( fcount * sizeof( struct struc_func ) ) ;
    d = ( struct struc_func* )malloc( dcount * sizeof( struct struc_func ) ) ;
    if( fread( f , sizeof( struct struc_func ) , fcount , fabc_dabc ) != fcount ) return ;
    if( fread( d , sizeof( struct struc_func ) , dcount , fabc_dabc ) != dcount ) return ;
  }
  fclose( fabc_dabc ) ;
#else // wrap to the computation
  actually_compute_fs_and_ds( ) ;
#endif 
  return ;
}

// contraction with I.d^{abc}A^{a}B^{b}C^{c}
double complex
dabc_ABC( const GLU_complex A[ NCNC ] , 
	  const GLU_complex B[ NCNC ] , 
	  const GLU_complex C[ NCNC ] )
{
  GLU_complex a[ NCNC-1 ] , b[ NCNC-1 ] , c[ NCNC-1 ] ;
  lie_data( a , A ) ;
  lie_data( b , B ) ;
  lie_data( c , C ) ;
  double complex result = 0.0 ;
  size_t i ;
  for( i = 0 ; i < dcount ; i++ ) {
    const int a_idx = d[i].a , b_idx = d[i].b , c_idx = d[i].c ;
    result += ( d[i].val ) * a[a_idx] * b[b_idx] * c[c_idx] ;
  }
  return result ;
}

// free the f and d matrices
void
free_f_and_d( void )
{
  free( f ) ;
  free( d ) ;
  return ;
}

// free lambda
void
free_generators( void )
{
  size_t i ;
  for( i = 0 ; i < NCNC-1 ; i++ ) { free( lambda[i] ) ; }
  free( lambda ) ;
  return ;
}

// contraction with I.f^{abc}A^{a}B^{b}C^{c}
double complex
ifabc_ABC( const GLU_complex A[ NCNC ] , 
	   const GLU_complex B[ NCNC ] , 
	   const GLU_complex C[ NCNC ] )
{
  GLU_complex a[ NCNC-1 ] , b[ NCNC-1 ] , c[ NCNC-1 ] ;
  lie_data( a , A ) ;
  lie_data( b , B ) ;
  lie_data( c , C ) ;
  double complex result = 0.0 ;
  size_t i ;
  for( i = 0 ; i < fcount ; i++ ) {
    const size_t a_idx = f[i].a , b_idx = f[i].b , c_idx = f[i].c ;
    result += ( I * f[i].val ) * a[a_idx] * b[b_idx] * c[c_idx] ;
  }
  return result ;
}

// contraction with the combination A^a B^b C^c (I.f^{abc}+d^{abc}) == 4.Tr(ABC)
double complex
ifabc_dabc_ABC( const GLU_complex A[ NCNC ] , 
		const GLU_complex B[ NCNC ] , 
		const GLU_complex C[ NCNC ] ) 
{
  // lie elements
  GLU_complex a[ NCNC-1 ] , b[ NCNC-1 ] , c[ NCNC-1 ] ;
  lie_data( a , A ) ;
  lie_data( b , B ) ;
  lie_data( c , C ) ;
  double complex result = 0.0 ;
  size_t i ;
  for( i = 0 ; i < dcount ; i++ ) {
    const size_t a_idx = d[i].a , b_idx = d[i].b , c_idx = d[i].c ;
    result += ( d[i].val ) * a[a_idx] * b[b_idx] * c[c_idx] ;
  }
  for( i = 0 ; i < fcount ; i++ ) {
    const size_t a_idx = f[i].a , b_idx = f[i].b , c_idx = f[i].c ;
    result += ( I*f[i].val ) * a[a_idx] * b[b_idx] * c[c_idx] ;
  }
  return result ;
}

// initialise generators of fundamental SU(NC)
void
init_generators( void )
{
  lambda = ( GLU_complex** )malloc( ( NCNC-1 ) * sizeof( GLU_complex* ) ) ;
  size_t i , j ;
  for( i = 0 ; i < NCNC-1 ; i++ ) {
    lambda[ i ] = ( GLU_complex* )malloc( NCNC * sizeof( GLU_complex ) ) ;
    for( j = 0 ; j < NCNC ; j++ ) { lambda[i][j] = 0.0 ; } // initialise to 0
  }  
  // the idea is to compute the 2x2, and then embed in the 3x3 and so on
  size_t subblock = 2 , idx = 0 ;
  while( subblock <= NC ) {
    // fix the column to be the last element in the subblock
    const size_t k = subblock-1 ;
    for( j = 0 ; j < k ; j++ ) {
      lambda[idx][ mat_idx(j,k) ] = 0.5 ;
      lambda[idx][ mat_idx(k,j) ] = 0.5 ;
      idx ++ ;
      lambda[idx][ mat_idx(j,k) ] = I * 0.5 ;
      lambda[idx][ mat_idx(k,j) ] = -I * 0.5 ;
      idx ++ ;
    } // fill in the final diagonal for this sub-block
    const register GLU_real fact = sqrt( 1.0 / (GLU_real)( 2.0 * (subblock) * ( subblock-1 ) ) ) ; 
    for( i = 0 ; i < subblock-1 ; i++ ) {
      lambda[idx][ mat_idx(i,i) ] = fact ;
    }
    lambda[idx][ mat_idx(i,i) ] = fact * ( -i ) ;
    idx++ ;
    // increment the size of the block matrix we are in
    subblock ++ ;
  }
  return ;
}
