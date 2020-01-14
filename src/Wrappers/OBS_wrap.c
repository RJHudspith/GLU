/*
    Copyright 2013-2018 Renwick James Hudspith

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
#include "invert.h"

// some extras we don't usually print out
//#define GLU_LOG_SPEED
//#define GLU_PRINT_LOOPS

// log types we are looking at
typedef enum { EXACT_FAST , EXACT_SLOW ,
	       EXACT_SLOWER , LINEAR } log_type ;

#if !(defined NO_LINK_CHECK)
// checks for any problems with non-SU( NC ) matrices
static void
check_links( const struct site *__restrict lat )
{
  const GLU_bool FUCKED = GLU_FALSE ;
  size_t i , NBADLINKS = 0 ;
  fprintf( stdout , "\n[UNITARY] Checking for Unitary link matrices ...\n" ) ;
  #pragma omp parallel for private(i) reduction(+:NBADLINKS)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu ; 
    for( mu = 0 ; mu < ND ; mu++ ) {
      if( is_unitary( lat[i].O[mu] ) == FUCKED ) {
	fprintf( stderr , "BAD LINK site :: %zu mu:: %zu \n", i , mu ) ;
	fprintf( stderr , "Determinant :: " ) ;
	printcomplex( det( lat[i].O[mu] ) ) ;
	NBADLINKS = NBADLINKS + (int)1 ;
      }
    }
  }
  if( NBADLINKS > 0 ) {
    fprintf( stderr , "[UNITARY] BADLINKS :: %zu out of %zu \n\n", 
	     NBADLINKS , LVOLUME * ND ) ;
  } else {
    fprintf( stdout , "[UNITARY] BADLINKS :: %zu out of %zu \n\n", 
	     NBADLINKS , LVOLUME * ND ) ;
  }
  return ;
}
#endif

#ifdef GLU_PRINT_LOOPS
/// 1x1 wilson loop check
static double
complete_plaquette( const GLU_complex *__restrict a ,
		    const GLU_complex *__restrict b ,
		    const GLU_complex *__restrict c , 
		    const GLU_complex *__restrict d )
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
  fprintf( stdout , "[GFACC] Looking at Landau Gauge-Fixing Information \n" ) ; 
  fprintf( stdout , "[GFACC LANDAU]      LIN :: %1.6e \n" , 
	   theta_test_lin( lat , &max , ND ) ) ;  
  fprintf( stdout , "[GFACC MAX LANDAU]  LIN :: %1.6e \n" , max ) ; 
#if NC < 4
  fprintf( stdout , "[GFACC LANDAU]      LOG :: %1.6e \n" , 
	   theta_test_log( lat , &max , ND ) ) ;  
  fprintf( stdout , "[GFACC MAX LANDAU]  LOG :: %1.6e \n" , max ) ; 
  double lin , log ;
  const_time( lat , &lin , &log ) ; 
  fprintf( stdout , "[GFACC] Temporal constance || LIN %e || Log %e " , 
	   lin , log ) ;
#endif
  // look at the coulomb stuff
  fprintf( stdout , "\n[GFACC] Looking at Coulomb "
	   "Gauge-Fixing Information \n" ) ; 
  fprintf( stdout , "[GFACC COULOMB]     LIN :: %1.6e \n" , 
	  theta_test_lin( lat , &max , ND - 1 ) ) ;  
  fprintf( stdout , "[GFACC MAX COULOMB] LIN :: %1.6e \n" , max ) ; 
#if NC < 4
  fprintf( stdout , "[GFACC COULOMB]     LOG :: %1.6e \n" , 
	  theta_test_log( lat , &max , ND - 1 ) ) ;  
  fprintf( stdout , "[GFACC MAX COULOMB] LOG :: %1.6e \n" , max ) ; 
#endif
  fprintf( stdout , "\n") ;
  return;
}

#ifdef GLU_PRINT_LOOPS
// output the invariance ...
static void
print_invariant_loops( const struct site *__restrict lat )
{
  size_t i ; 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu , nu , t , s ;
    for( mu = 1 ; mu < ND ; mu++ ) {
      t = lat[i].neighbor[mu] ; 
      for( nu = 0 ; nu < mu ; nu++ ) {
	s = lat[i].neighbor[nu] ;
	const double face = complete_plaquette( lat[ i ].O[mu] , 
						lat[ t ].O[nu] , 
						lat[ s ].O[mu] , 
						lat[ i ].O[nu] ) ; 
	fprintf( stdout , " %d %d %d :: %1.15f \n" , i , mu , nu , face ) ; 
      }
    }
  }
  return ;
}
#endif

#ifdef GLU_LOG_SPEED
// test the log speed
static void
test_log_speed( const struct site *__restrict lat , 
		const log_type log )
{
  GLU_complex temp2[ NCNC ] ;
  size_t maxstress = 250 , stress , elements ;
  double times[ maxstress ] ;
  double sum = 0.0 ;
  size_t mu , i ;
  // stress loops
  for( stress = 0 ; stress < maxstress ; stress++ ) {
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
    fprintf( stdout , "[LOG ERROR] %e \n", sum / ( ND * NCNC ) ) ;
  }
  double ave = 0.0 , err = 0.0 ;
  average( &ave , &err , times , stress ) ;
  // compute the average of times ...
  fprintf( stdout , "RESULT :: %e %e \n" , ave , err ) ;
  return ;
}
#endif

// wrapper for the default behaviour
void
gauge( const struct site *__restrict lat )
{
  start_timer( ) ;

  fprintf( stdout , "\n" ) ;
  double splaq , tplaq , plaq = all_plaquettes( lat , &splaq , &tplaq ) ;
  fprintf( stdout , "[PLAQS SU(%d)]          :: %1.15f \n" , NC , plaq ) ;
  fprintf( stdout , "[PLAQS SU(%d) spatial]  :: %1.15f \n" , NC , splaq ) ;
  fprintf( stdout , "[PLAQS SU(%d) temporal] :: %1.15f \n" , NC , tplaq ) ;
  fprintf( stdout , "\n" ) ;
  double slink , tlink , link = all_links( lat , &slink , &tlink ) ;
  fprintf( stdout , "[LINKS SU(%d)]          :: %1.15f \n" , NC , link ) ;
  fprintf( stdout , "[LINKS SU(%d) spatial]  :: %1.15f \n" , NC , slink ) ;
  fprintf( stdout , "[LINKS SU(%d) temporal] :: %1.15f \n" , NC , tlink ) ;
  fprintf( stdout , "\n") ;

  // gauge invariance checks and what have you
#if !(defined NO_LINK_CHECK)

  printf( "WHAT!!!!!!!!!!!!!!!\n" ) ;
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // compute correct denominator
    size_t denom = 1 , nu ;
    for( nu = 0 ; nu < ND ; nu++ ) {
      denom *= ( nu != mu ) ? Latt.dims[nu] : 1 ;
    }
    // average
    const double DENOM = 1.0 / ( ( double )NC * denom ) ;
    double complex R = poly( lat , mu ) ;
    fprintf( stdout , "[POLY SU(%d) DIR_%zu] ( %e , %e ) , MODULUS :: %e \n" ,
	    NC , mu , creal( R ) * DENOM , cimag( R ) * DENOM , 
	    cabs( R ) * DENOM ) ;
  }
  check_links( lat ) ;
#endif

  // check the gauge fixing accuracy for the LINEAR and LOG definition
  gf_check( lat ) ;

  print_time( ) ;

  #ifdef GLU_LOG_SPEED
    // have a look at various log speeds and what have you
    test_log_speed( lat , EXACT_SLOW ) ;
  #elif defined GLU_PRINT_LOOPS
    // if we want to look at the 1x1 plaquette invariance
    print_invariant_loops( lat ) ;
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
  #undef GLU_TEST_PRODS
#endif
