/*
    Copyright 2013 Renwick James Hudspith

    This file (GF_wrap.c) is part of GLU.

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
   @file GF_wrap.c
   @brief rappers for the gauge fixing aspects of the code.
 */

// is special
#include "Mainfile.h"

#include "chklat_stuff.h" // wanted for the skip_hdr function ...
#include "Coulomb.h"      // Coulomb fixing
#include "gtrans.h"       // gauge transformations
#include "GLU_timer.h"    // for the timer
#include "GLU_memcheck.h" // to tell us if we have the memory
#include "Landau.h"       // Landau fixing
#include "MAG.h"          // Maximal Axial Gauge
#include "Or.h"     // Landau fixing
#include "plaqs_links.h"  // for the plaquettes and links
#include "read_config.h"  // read the config back in
#include "read_headers.h" // for rereading the header
#include "SM_wrap.h"      // smeared preconditioning

#ifdef LUXURY_GAUGE
#include "geometry.h"
#endif

// coulomb wrapper
static int 
GF_wrap_coulomb( struct site *__restrict lat ,
		 const struct gf_info GFINFO )
{
  start_timer( ) ;

  // if we don't have the memory we leave
  if( have_memory_Cgf( ) == GLU_FAILURE ) { return GLU_FAILURE ; }

  if( GFINFO.improve == MAG_IMPROVE ) {
    #if !(defined deriv_full || defined deriv_fulln )
    /// ewww! have to both malloc and free gauge here -> Memory dangerous?
    GLU_complex **gauge = malloc( LVOLUME * sizeof( GLU_complex* ) ) ; 
    int i ;
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      GLU_malloc( (void**)&gauge[i] , 16 , NCNC * sizeof( GLU_complex ) ) ;
    }
    // mag improvement
    mag( lat , gauge ) ; 
    //axial_gauge( lat , gauge , ND-1 ) ;
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      free( gauge[i] ) ;   
    }
    free( gauge ) ;  
    #else
    fprintf( stderr , "MAG not compatible with exact-log fixing \n" ) ;
    #endif
  }

#ifdef OVERRELAXED_GF
  // we could check iters if we wanted, actually we do want to
  double acc ;
  OrCoulomb( lat , &acc , GFINFO.max_iters ,
	     GFINFO.accuracy , Latt.gf_alpha ) ;
#else
  // we could check iters if we wanted, actually we do want to
  Coulomb( lat , GFINFO.accuracy , GFINFO.max_iters ) ; 
#endif

  // fixes the average temporal links to 1 as best it can
  // using an analogue of the temporal gauge ...
  if( GFINFO.improve == RESIDUAL_IMPROVE ) {
    fprintf( stdout , "\nResidual gauge fixing step \n" ) ;
    residual_fix( lat ) ;
    fprintf( stdout , "[TLINK] %1.15f [SLINK] %1.15f "
	     "[LINK] %1.15f [PLAQ] %1.15f \n" , 
	     t_links(lat) , s_links(lat) , links(lat) , av_plaquette(lat) ) ;
  }
  print_time( ) ;

  return GLU_SUCCESS ;
}

// landau wrapper
static int
GF_wrap_landau( struct site *__restrict lat ,
		const char *__restrict infile ,
		const struct gf_info GFINFO ,
		const GF_improvements improvement )
{
  int iters = 0 ;
  start_timer( ) ;
#ifdef OVERRELAXED_GF
  // have to alloc a temporary gauge field for this
  if( GFINFO.improve == MAG_IMPROVE ) {
    GLU_complex **gauge = NULL ;
    if( GLU_malloc( (void**)&gauge , 16 , LVOLUME * sizeof( GLU_complex* ) ) != 0 ) {
      fprintf( stderr , "[GF] GF_wrap_landau (overrelaxed-MAG) failed "
	       "to allocate temporary gauge\n" ) ;
      return GLU_FAILURE ;
    }
    size_t i ;
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      GLU_malloc( (void**)&gauge[i] , 16 , NCNC * sizeof( GLU_complex ) ) ;
      identity( gauge[i] ) ;
    }
    mag( lat , gauge ) ;
    // free the gauge transformation matrices
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      free( gauge[i] ) ;   
    }
    free( gauge ) ;
  }
  // Overrelaxed gauge fixing routine
  double acc ;
  iters = OrLandau( lat , &acc , GFINFO.max_iters ,
		    GFINFO.accuracy , Latt.gf_alpha ) ;
#else
  GLU_complex **gauge = NULL ;
  if( GLU_malloc( (void**)&gauge , 16 , LVOLUME * sizeof( GLU_complex* ) ) != 0 ) {
    fprintf( stderr , "[GF] GF_wrap_landau failed to allocate "
	     "temporary gauge\n" ) ;
    return GLU_FAILURE ;
  }
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    gauge[i] = ( GLU_complex* )malloc( NCNC * sizeof( GLU_complex ) ) ; 
    identity( gauge[i] ) ;
  }

  // start the timer unless we are smeared-preconditioning
  if( improvement != SMPREC_IMPROVE ) { start_timer( ) ; }

  // MAG - preconditioning..
  if( GFINFO.improve == MAG_IMPROVE ) { 
    #if !(defined deriv_full || defined deriv_fulln )
    mag( lat , gauge ) ; 
    #else
    fprintf( stderr , "MAG not compatible with exact-log fixing \n" ) ;
    #endif
  }  

  // the memory cheap one wasn't much of a saving so we just use the fast
  iters = Landau( lat , gauge , 
		  GFINFO.accuracy , 
		  GFINFO.max_iters , 
		  infile , improvement ) ; 
  // free the gauge transformation matrices
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    free( gauge[i] ) ;   
  }
  free( gauge ) ;
#endif
  print_time( ) ;
  return iters ;
}

// print out the details
static void
print_fixing_info( const struct gf_info GFINFO ,
		   const struct sm_info SMINFO )
{
#ifdef OVERRELAXED_GF
  fprintf( stdout , "\n[GF] Using the Over-Relaxation routines\n[GF] " ) ;
#else
#ifdef GLU_GFIX_SD
  #ifdef HAVE_FFTW3_H
  fprintf( stdout , "\n[GF] Using the Fourier accelerated steepest "
	   "descent (FASD) routines\n[GF] " ) ;
  #else
  fprintf( stdout , "\n[GF] Using the steepest descent (SD) routines\n[GF] " ) ;
  #endif
#else // default is the CG
  #ifdef HAVE_FFTW3_H
  fprintf( stdout , "\n[GF] Using the Fourier accelerated non-linear"
	   " conjugate gradient (FACG) routines\n[GF] " ) ;
  #else
  fprintf( stdout , "\n[GF] Using the non-linear"
	   " conjugate gradient (CG) routines\n[GF] " ) ;
  #endif
#endif
#endif
  switch( GFINFO.type ) {
  case GLU_LANDAU_FIX :
    fprintf( stdout , "Fixing to Landau gauge, within an accuracy of %g \n" ,
	     GFINFO.accuracy ) ; 
    break ;
  case GLU_COULOMB_FIX :
    fprintf( stdout , "Fixing to Coulomb gauge, within an accuracy of %g \n" ,
	     GFINFO.accuracy ) ; 
    break ;
  default : break ;// should never get here
  }
  fprintf( stdout , "[GF] Tuning parameter alpha :: %f \n" , Latt.gf_alpha ) ;
  fprintf( stdout , "[GF] Performing AT MOST %zu iterations"
	   "before randomly restarting ... \n" , GFINFO.max_iters ) ;
  fprintf( stdout , "[GF] Allowing for %d restarts before complaint ... \n" , 
	   GF_GLU_FAILURES ) ; 
  // derivative routines available
  fprintf( stdout , "[GF] " ) ;
  #ifdef deriv_lin
  fprintf( stdout , "Using the Hermitian projection definition"
	   "of the gluon fields.\n" ) ; 
  #elif defined deriv_linn
  fprintf( stdout , "Using the Hermitian projection definition"
	   "of the gluon fields,\n"
	   "and the next nearest neighbor derivative.\n" ) ; 
  fprintf( stdout , "[GF] Nearest neighbour terms nn1 :: %f nn2 :: %f \n" , 
	   nn1 , nn2 ) ;
  #elif defined deriv_full
  fprintf( stdout , "Using the exact Log definition of the gluon fields.\n" ) ; 
  #elif defined deriv_fulln
  fprintf( stdout , "Using the exact Log definition of the gluon fields,\n"
	   "and the next nearest neighbor derivative.\n" ) ;
  fprintf( stdout , "[GF] Nearest neighbour terms nn1 :: %f nn2 :: %f \n" , 
	   nn1 , nn2 ) ;
  #elif defined deriv_fullnn
  fprintf( stdout , "Using the exact Log definition of the gluon fields,\n"
	   "and the next to next nearest neighbor derivative improvement.\n" ) ;
  #endif
  // exponentiation approximation routines
  fprintf( stdout , "[GF] " ) ;
  #if defined exp_approx
  fprintf( stdout , "Approximate O(a) exponential expansion,"
	   "and reunitarisation in Gauge Fixing.\n" ) ; 
  #elif defined exp_a2_approx
  fprintf( stdout , "Approximate O(a^2) exponential expansion,"
	   "and reunitarisation in Gauge Fixing.\n" ) ; 
  #elif defined exp_exact
  fprintf( stdout , "Exact exponentiation in the Gauge Fixing being used.\n" ) ;
  #endif
// tell us which log-method we are using
#if ( defined deriv_full ) || ( defined deriv_fulln )
  fprintf( stdout , "[GF] Using Vandermonde logarithmic def warm-up \n" ) ;
#endif
#ifdef LUXURY_GAUGE
  fprintf( stdout , "[GF] " ) ;
  #ifdef WORST_COPY
  fprintf( stdout , "Finding the worst from %d random copies \n" , 
	   LUXURY_GAUGE ) ;
  #else
  fprintf( stdout , "Finding the best from %d random copies \n" , 
	   LUXURY_GAUGE ) ;
  #endif
#endif
  fprintf( stdout , "[GF] " ) ;
  switch( GFINFO.improve ) {
  case MAG_IMPROVE :
    fprintf( stdout , "MAG-Preconditioning Improvement.\n" ) ; 
    break ;
  case SMPREC_IMPROVE :
    if( GFINFO.type == GLU_COULOMB_FIX )
      fprintf( stderr , "Smearing-Preconditioning not implemented.\n" ) ;
    else
      fprintf( stdout , "Smearing-Preconditioning Improvement.\n" ) ; 
    break ;
  case RESIDUAL_IMPROVE :
    if( GFINFO.type != GLU_LANDAU_FIX ) {
      fprintf( stdout , "Residual gauge fixing improvement, "
	       "after Coulomb fixing.\n" ) ;
    }
    break ;
  default :
    fprintf( stdout , "No improvement.\n" ) ;
    break ;
  }
  fprintf( stdout , "\n" ) ;
  return ;
}

// smeared preconditioned step
static int
smeared_prec_step( lat , gauge , SMINFO , HEAD_DATA , infile , Local_maxiters , Local_accuracy , improvement )
     struct site *__restrict lat ;
     GLU_complex *__restrict *__restrict gauge ;
     const struct sm_info SMINFO ; 
     struct head_data HEAD_DATA ;
     const char *__restrict infile ;
     const int Local_maxiters ;
     const double Local_accuracy ;
     const int improvement ;
{
  if( SM_wrap_struct( lat , SMINFO ) == GLU_FAILURE ) { return GLU_FAILURE ; } 

  // we set the alpha to be smaller as with smooth fields the tuning graph is different!
  Latt.gf_alpha *= 7./8. ; 
  //Landau_4smear( lat , gauge , Local_accuracy , Local_maxiters ) ; 
  Landau( lat , gauge , Local_accuracy , Local_maxiters , infile , SMPREC_IMPROVE ) ; 
  // set it back to normal
  Latt.gf_alpha *= 8./7. ; 

  // this should be wrapped into grab file no?
  struct head_data temp_head ;
  FILE *config = fopen( infile , "rb" ) ;
  if( read_header( config , &temp_head , GLU_FALSE ) == GLU_FAILURE ) return GLU_FAILURE ; 
  if( get_config_SUNC( config , lat , temp_head ) == GLU_FAILURE ) return GLU_FAILURE ;
  gtransform( lat , (const GLU_complex **)gauge ) ; 

  return GLU_SUCCESS ;
}

#ifdef LUXURY_GAUGE

////////// LUXURY GAUGE SMPREC VERSION /////////////////
/// leaves the best copy in "lat"
static int
GF_wrap_smprec_luxury( struct site *__restrict lat ,
		       const char *__restrict infile ,
		       const struct gf_info GFINFO ,
		       const struct sm_info SMINFO ,
		       const struct head_data HEAD_DATA )
{
  // need to allocate a temporary lattice with this one ...
  struct site *lat_cpy = NULL ;
  GLU_complex **gauge = NULL ;
  
  if( GLU_malloc( (void**)&gauge   , 16 ,  LVOLUME * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&lat_cpy , 16 ,  LVOLUME * sizeof( struct site )  ) != 0 ) {
    fprintf( stderr , "[GF] GF_wrap_smprec_luxury failed "
	     "to allocate temporary space\n" ) ;
    return GLU_FAILURE ;
  }
  init_navig( lat_cpy ) ;

  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    gauge[i] = ( GLU_complex* )malloc( NCNC * sizeof( GLU_complex ) ) ; 
  }

  // set up a local number of iterations and a local smearing accuracy why not?
  const int Local_maxiters = 1800 ;
  const double Local_accuracy = 1E-8 ;

  double MIN_GFUNC = 1.0 ;

  // loop over the number of gribov copies ...
  int copies , cpy_idx = 0 ;
  for( copies = 0 ; copies < LUXURY_GAUGE ; copies++ ) {
    // read in the lattice an perform a random transform
    if( grab_file( lat_cpy , gauge , infile ) == GLU_FAILURE ) { 
      return GLU_FAILURE ; }

    smeared_prec_step( lat_cpy , gauge , SMINFO , HEAD_DATA , 
		       infile , Local_maxiters , Local_accuracy , 
		       GFINFO.improve ) ;

    // set gauge to the identity
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) { identity( gauge[i] ) ; }

    // this can be done roughly too
    const int iters = Landau( lat_cpy , gauge , 
			      Local_accuracy , Local_maxiters , 
			      infile , SMPREC_IMPROVE ) ; 

    // compute the gauge functional
    const double GFUNC = links( lat ) ;
    fprintf( stdout , "\n   [COPY] %d [FUNCTIONAL] %1.15f [ITER] %d " , 
	     copies , GFUNC , iters ) ; 
    if( GFUNC < MIN_GFUNC && iters != Local_maxiters ) {
      fprintf( stdout , " -> Copy accepted \n" ) ;
      // copy them over
      MIN_GFUNC = GFUNC ;
      cpy_idx = copies ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	memcpy( &lat[i] , &lat_cpy[i] , sizeof( struct site ) ) ;
      }
    } else { fprintf( stdout , " -> Copy rejected \n\n" ) ; }
  print_time( ) ;
  fprintf( stdout , "\n" ) ;
  }

  // completion run ...
  Landau( lat , gauge , 
	  GFINFO.accuracy , GFINFO.max_iters , 
	  infile , NO_IMPROVE ) ; 

  // the output is a little confusing so I print out the important details here
  fprintf( stdout , "\n[GAUGE_COPY SELECTED] %d "
	   "[FUNCTIONAL] %1.15f [PLAQUETTE] %1.15f \n" , 
	   cpy_idx , MIN_GFUNC , av_plaquette(lat) ) ;
  GLU_real tr ;
  const double link = indivlinks( lat , &tr ) ;
  fprintf( stdout , "[LINK] %1.15f [MAX] %1.15f \n" , link , tr/(GLU_real)NC ) ;

#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    free( gauge[i] ) ; 
  }
  free( gauge ) ;
  free( lat_cpy ) ;

  return GLU_SUCCESS ;
}

#else

// smeared-preconditioned step
static int
GF_wrap_smprec( struct site *__restrict lat ,
		const char *__restrict infile ,
		const struct gf_info GFINFO ,
		const struct sm_info SMINFO ,
		const struct head_data HEAD_DATA ,
		int recurses )
{
  GLU_complex **gauge = NULL ;
  if( GLU_malloc( (void**)&gauge , 16 , LVOLUME*sizeof( GLU_complex ) ) != 0 ) {
    fprintf( stderr , "[GF] GF_wrap_smprec failed to "
	     "allocate temporary fields \n" ) ;
    return GLU_FAILURE ;
  }
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    gauge[i] = ( GLU_complex* )malloc( NCNC * sizeof( GLU_complex ) ) ; 
  }
  // set up a local number of iterations and a local smearing accuracy why not?
  const int Local_maxiters = 1500 ;
  const double Local_accuracy = 1E-8 ;

  smeared_prec_step( lat , gauge , SMINFO , HEAD_DATA , 
		     infile , Local_maxiters , Local_accuracy ,
		     SMPREC_IMPROVE ) ;

  // if this fails we reread and call the whole thing again
  GLU_bool restart = GLU_FALSE ;
  if( GF_wrap_landau( lat , infile , GFINFO , NO_IMPROVE ) == 
      GFINFO.max_iters && 
      ( recurses < GF_GLU_FAILURES ) ) {
    if( grab_file( lat , gauge , infile ) == GLU_FAILURE ) {
      fprintf( stderr , "[GF] Error in (re)reading the file ..."
	       "Leaving in a state of disarray\n" ) ;
    } else {
      restart = GLU_TRUE ;
    }
  } if( recurses == GF_GLU_FAILURES ) {
    fprintf( stderr , "\n[GF] We have failed enough! %d Strongly consider"
	     "increasing the tuning parameter"
	    " and/or increasing the number of gauge-fixing iterations.\n" , 
	     recurses ) ;
  }
  print_time( ) ;

  // free the gauge transformation matrices
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    free( gauge[i] ) ; 
  }
  free( gauge ) ;
  
  // here is the call
  if( restart == GLU_TRUE ) {
    fprintf( stdout , "\n\n ************ GLU_FAILURE %d **************** \n\n" ,
	     recurses ) ;
    GF_wrap_smprec( lat , infile , GFINFO , SMINFO , HEAD_DATA , ++recurses ) ;
  }

  return GLU_SUCCESS ;
}
#endif

// have a wrapper for the above functions here ...
int 
GF_wrap( const char *infile , 
	 struct site *__restrict lat , 
	 const struct gf_info GFINFO , 
	 const struct sm_info SMINFO ,
	 const struct head_data HEAD_DATA )
{
  if( GFINFO.type == DEFAULT_NOFIX ) return GLU_SUCCESS ;
  print_fixing_info( GFINFO , SMINFO ) ;
  // setup for Landau - Coulomb
  if( GFINFO.type == GLU_LANDAU_FIX ) {
    if( GFINFO.improve == SMPREC_IMPROVE ) {
      #ifdef LUXURY_GAUGE
      return GF_wrap_smprec_luxury( lat , infile , GFINFO , 
				    SMINFO , HEAD_DATA ) ;
      #else
      int recurses = 0 ;
      return GF_wrap_smprec( lat , infile , GFINFO , SMINFO , 
			     HEAD_DATA , recurses ) ;
      #endif 
    } else {
      return GF_wrap_landau( lat , infile , GFINFO , GFINFO.improve ) ;
    }
  } else if( GFINFO.type == GLU_COULOMB_FIX ) {
    return GF_wrap_coulomb( lat , GFINFO ) ;
  }
  return GLU_SUCCESS ;
}
