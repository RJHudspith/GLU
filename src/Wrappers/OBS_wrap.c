/*
    Copyright 2013 Renwick James Hudspith

    This file (OBS_wrap.c) is part of GLU.

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
   @file OBS_wrap.c
   @brief wraps many observable calculations into one

   calls av_plaquette() is_unitary() poly() s_plaq() t_plaq() links() s_links() t_links()
 */

#include "Mainfile.h"

#include "clover.h"
#include "gftests.h"
#include "GLU_timer.h"
#include "plaqs_links.h"
#include "POLY.h"

// some extras we don't usually print out
//#define GLU_LOG_SPEED
//#define GLU_PRINT_LOOPS
//#define GLU_TEST_MULS
//#define GLU_TEST_PRODS

// log types we are looking at
typedef enum { EXACT_FAST , EXACT_SLOW ,
	       EXACT_SLOWER , LINEAR } log_type ;

// checks for any problems with non-SU( NC ) matrices
static void
check_links( lat )
     const struct site *__restrict lat ;
{
  const GLU_bool FUCKED = GLU_TRUE ;
  int i , NBADLINKS = 0 ;
  printf( "\n[UNITARY] Checking for Unitarity of link matrices ...\n" ) ;
#pragma omp parallel for private(i) reduction(+:NBADLINKS)
  for( i = 0 ; i < Latt.Volume ; i++ ) {
    int mu ; 
    for( mu = 0 ; mu < ND ; mu++ ) {
      if( unlikely( is_unitary( lat[i].O[mu] ) == FUCKED ) ) {
	printf( "BAD LINK site :: %d mu:: %d \n", i , mu ) ;
	printf( "Determinant :: " ) ;
	printcomplex( det( lat[i].O[mu] ) ) ;
	NBADLINKS = NBADLINKS + (int)1 ;
      }
    }
  }
  printf( "[UNITARY] BADLINKS :: %d out of %d \n\n", 
	  NBADLINKS , LVOLUME * ND ) ;
  return ;
}

#ifdef GLU_PRINT_LOOPS
/// 1x1 wilson loop check
static double
complete_plaquette( a , b , c , d )
     const GLU_complex *__restrict a , *__restrict b ;
     const GLU_complex *__restrict c , *__restrict d ;
{
  GLU_complex temp1[ NCNC ] , temp2[ NCNC ] ; 
  multab( temp1 , a , b ) ; 
  multab_dag( temp2 , temp1 , c ) ; 
  multab_dag( temp1 , temp2 , d ) ; 
  double tr ;
  speed_trace_Re( &tr , temp1 ) ;
  return tr ;
}
#endif

// check the gauge fixing of the config
static void
gf_check( const struct site *__restrict lat )
{
  GLU_real max ; 
  // landau information
  printf( "[GFACC] Looking at Landau Gauge-Fixing Information \n" ) ; 
  printf( "[GFACC LANDAU]      LIN :: %1.6e \n" , 
	  theta_test_lin( lat , &max , ND ) ) ;  
  printf( "[GFACC MAX LANDAU]  LIN :: %1.6e \n" , max ) ; 
#if NC < 4
  printf( "[GFACC LANDAU]      LOG :: %1.6e \n" , 
	  theta_test_log( lat , &max , ND ) ) ;  
  printf( "[GFACC MAX LANDAU]  LOG :: %1.6e \n" , max ) ; 
  double lin , log ;
  const_time( lat , &lin , &log ) ; 
  printf( "[GFACC] Temporal constance || LIN %e || Log %e " , lin , log ) ;
#endif
  // look at the coulomb stuff
  printf( "\n[GFACC] Looking at Coulomb Gauge-Fixing Information \n" ) ; 
  printf( "[GFACC COULOMB]     LIN :: %1.6e \n" , 
	  theta_test_lin( lat , &max , ND - 1 ) ) ;  
  printf( "[GFACC MAX COULOMB] LIN :: %1.6e \n" , max ) ; 
#if NC < 4
  printf( "[GFACC COULOMB]     LOG :: %1.6e \n" , 
	  theta_test_log( lat , &max , ND - 1 ) ) ;  
  printf( "[GFACC MAX COULOMB] LOG :: %1.6e \n" , max ) ; 
#endif
  printf("\n") ;
  return;
}

#ifdef GLU_PRINT_LOOPS
// output the invariance ...
static void
print_invariant_loops( lat )
     const struct site *__restrict lat ;
{
  int i ; 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    int mu ;
    for( mu = 1 ; mu < ND ; mu++ ) {
      register const int t = lat[i].neighbor[mu] ; 
      int nu ;
      for( nu = 0 ; nu < mu ; nu++ ) {
	register const int s = lat[i].neighbor[nu] ;
	const double face = complete_plaquette( lat[ i ].O[mu] , 
						lat[ t ].O[nu] , 
						lat[ s ].O[mu] , 
						lat[ i ].O[nu] ) ; 
	printf( " %d %d %d :: %1.15f \n" , i , mu , nu , face ) ; 
      }
    }
  }
  return ;
}
#endif

#ifdef GLU_LOG_SPEED

// stabler average computation ...
static void
average( double *ave , double *err , 
	 const double *data , const int N )
{
  int i ;
  *err = *ave = 0.0 ;
  for( i = 0 ; i < N ; i++ ) {
    const double delta = data[i] - *ave ;
    *ave += delta / (double)i ;
    *err += delta * ( data[i] - *ave ) ;
  }
  *err /= ( N - 1 ) ;
}

// test the log speed
static void
test_log_speed( const struct site *__restrict lat , 
		const log_type log )
{
  GLU_complex temp2[ NCNC ] ;
  int maxstress = 250 , stress , elements ;
  double times[ maxstress ] ;
  double sum = 0.0 ;
  // stress loops
  for( stress = 0 ; stress < maxstress ; stress++ ) {
    int mu , i ;
    #ifdef TIME_GF
    start_timer() ;
    #endif
    for( i = 0 ; i < LVOLUME ; i++ ) {
      GLU_complex temp[ NCNC ] ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	switch( log ) {
	case EXACT_FAST : exact_log_fast( temp , lat[i].O[mu] ) ; break ;
	case EXACT_SLOW : exact_log_slow( temp , lat[i].O[mu] ) ; break ;
	case EXACT_SLOWER : Hermitian_proj( temp , lat[i].O[mu] ) ; break ;
	case LINEAR : get_iQ( temp , lat[i].O[mu] ) ; break ;
	}
	// exponentiation time is constant between them, ok to include ?
	exponentiate( temp2 , temp ) ;
	for( elements = 0 ; elements < NCNC ; elements++ ) {	    
	  sum += cabs( temp2[ elements ] -
		       lat[i].O[mu][ elements ] ) ;
	}
	// end of the line
      }
    }
    times[ stress ] = print_time() ;
    printf("[LOG ERROR] %e \n", sum / ( ND * NCNC ) ) ;
  }
  double ave = 0.0 , err = 0.0 ;
  average( &ave , &err , times , stress ) ;
  // compute the average of times ...
  printf( "RESULT :: %e %e \n" , ave , err ) ;
  return ;
}
#endif

// utility routine for testing matrix multiplies, just in case we change
// the matrix multiply routines
#ifdef GLU_TEST_MULS
static void
test_muls( const GLU_complex a[ NCNC ] ,
	   const GLU_complex b[ NCNC ] )
{
  GLU_complex res[ NCNC ] ;
  printf( "multab \n" ) ;
  multab_suNC( res , a , b ) ;
  write_matrix( res ) ;
  printf( "multabdag \n" ) ;
  multabdag_suNC( res , a , b ) ;
  write_matrix( res ) ;
  printf( "multabdagdag \n" ) ;
  multab_dagdag_suNC( res , a , b ) ;
  write_matrix( res ) ;
  printf( "multab_dag \n" ) ;
  multab_dag_suNC( res , a , b ) ;
  write_matrix( res ) ;
}
#endif

#ifdef GLU_TEST_PRODS
static void
test_trace_prods( const struct site *__restrict lat )
{
  printf( "[OBS] Trace products A.B \n" ) ;
  GLU_complex A[ NCNC ] , B[ NCNC ] ;
  Hermitian_proj( A , lat[1].O[0] ) ;
  Hermitian_proj( B , lat[1].O[1] ) ;

  GLU_complex temp[ NCNC ] ;
  multab( temp , A , B ) ;
  printcomplex( trace( temp ) ) ;
  
  double tr ;
  trace_ab_herm( &tr , A , B ) ;
  printf( " ( %1.7e  ,  %1.7e )\n" , tr , 0.0 ) ;

  GLU_complex Ctr ;
  trace_ab( &Ctr , A , B ) ;
  printcomplex( Ctr ) ;

  trace_ab_dag( &Ctr , A , B ) ;
  printcomplex( Ctr ) ;

  GLU_complex a[ HERMSIZE ] , b[ HERMSIZE ] ;
  Hermitian_proj_short( a , lat[1].O[0] ) ;
  Hermitian_proj_short( b , lat[1].O[1] ) ;

  trace_ab_herm_short( &tr , a , b ) ;
  printf( " ( %1.7e  ,  %1.7e )\n" , tr , 0.0 ) ;

  printf( "\n[OBS] Trace products A.A \n" ) ;

  multab( temp , A , A ) ;
  printcomplex( trace( temp ) ) ;
  
  trace_ab_herm_short( &tr , a , a ) ;
  printf( " ( %1.7e  ,  %1.7e )\n" , tr , 0.0 ) ;

  trace_prod_herm( &tr , A ) ;
  printf( " ( %1.7e  ,  %1.7e )\n" , tr , 0.0 ) ;

  printf( "\n[OBS] Trace products B.A.B \n" ) ;

  // product of three matrices can use the dag here
  // as it is the same through hermiticity
  multab_atomic_left( temp , B ) ;
  printcomplex( trace(temp) ) ;

  trace_abc( &Ctr , B , A , B ) ;
  printcomplex( Ctr ) ;

  trace_abc_dag( &Ctr , B , A , B ) ;
  printcomplex( Ctr ) ;

  trace_abc_dag_Re( &tr , B , A , B ) ;
  printf( " ( %1.7e  ,  %1.7e )\n" , tr , 0.0 ) ;
  return ;
}
#endif

// wrapper for the default behaviour
void
gauge( const struct site *__restrict lat )
{
  start_timer( ) ;

  printf("\n") ;
  double splaq , tplaq , plaq = all_plaquettes( lat , &splaq , &tplaq ) ;
  printf("[PLAQS SU(%d)]          :: %1.15f \n", NC , plaq ) ;
  printf("[PLAQS SU(%d) spatial]  :: %1.15f \n", NC , splaq ) ;
  printf("[PLAQS SU(%d) temporal] :: %1.15f \n", NC , tplaq ) ;
  printf("\n") ;
  double slink , tlink , link = all_links( lat , &slink , &tlink ) ;
  printf("[LINKS SU(%d)]          :: %1.15f \n", NC , link ) ;
  printf("[LINKS SU(%d) spatial]  :: %1.15f \n", NC , slink ) ;
  printf("[LINKS SU(%d) temporal] :: %1.15f \n", NC , tlink ) ;
  printf("\n") ;

  // gauge invariance checks and what have you
  int mu ;
  const double DENOM = 1.0 / ( ( double )NC * LCU ) ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    double complex R = poly( lat , mu ) ;
    printf( "[POLY SU(%d) DIR_%d] ( %e , %e ) , MODULUS :: %e \n" ,
	    NC , mu , creal( R ) * DENOM , cimag( R ) * DENOM , 
	    cabs( R ) * DENOM ) ;
  }
  check_links( lat ) ;

  // check the gauge fixing accuracy for the LINEAR and LOG definition
  gf_check( lat ) ;

  print_time( ) ;

  #ifdef GLU_LOG_SPEED
    // have a look at various log speeds and what have you
    test_log_speed( lat , EXACT_SLOW ) ;
  #elif defined GLU_PRINT_LOOPS
    // if we want to look at the 1x1 plaquette invariance
    print_invariant_loops( lat ) ;
  #elif defined GLU_TEST_MULS
    // check this against previous codes
    test_muls( lat[0].O[0] , lat[0].O[1] ) ;
  #elif defined GLU_TEST_PRODS
    test_trace_prods( lat ) ;
  #endif

  return ;
}

#ifdef GLU_LOG_SPEED
  #undef GLU_LOG_SPEED
#endif

#ifdef GLU_PRINT_LOOPS
  #undef GLU_PRINT_LOOPS
#endif

#ifdef GLU_TEST_MULS
  #undef GLU_TEST_MULS
#endif

#ifdef GLU_TEST_PRODS
  #undef GLU_TEST_MULS
#endif
