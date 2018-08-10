/*
    Copyright 2013-2018 Renwick James Hudspith

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
#include "Mainfile.h"

#include "chklat_stuff.h" // wanted for the skip_hdr function ...
#include "Coulomb.h"      // Coulomb fixing
#include "gtrans.h"       // gauge transformations
#include "GLU_timer.h"    // for the timer
#include "GLU_memcheck.h" // to tell us if we have the memory
#include "Landau.h"       // Landau fixing
#include "MAG.h"          // Maximal Axial Gauge
#include "Or.h"           // Landau fixing
#include "plaqs_links.h"  // for the plaquettes and links
#include "read_config.h"  // read the config back in
#include "read_headers.h" // for rereading the header
#include "SM_wrap.h"      // smeared preconditioning

// coulomb wrapper
static int 
GF_wrap_coulomb( struct site *lat ,
		 const struct gf_info GFINFO )
{
  start_timer( ) ;

  // if we don't have the memory we leave
  if( have_memory_Cgf( ) == GLU_FAILURE ) { return GLU_FAILURE ; }

#ifdef OVERRELAXED_GF
  // we could check iters if we wanted, actually we do want to
  double acc ;
  OrCoulomb( lat , &acc , GFINFO.max_iters ,
	     GFINFO.accuracy , Latt.gf_alpha ) ;
#else
  // we could check iters if we wanted, actually we do want to
  Coulomb( lat , GFINFO.accuracy , GFINFO.max_iters ) ; 
#endif

  print_time() ;

  return GLU_SUCCESS ;
}

// landau wrapper
static int
GF_wrap_landau( struct site *lat ,
		const char *infile ,
		const struct gf_info GFINFO )
{
  size_t iters = 0 ;
  start_timer( ) ;
#ifdef OVERRELAXED_GF
  // Overrelaxed gauge fixing routine
  double acc ;
  iters = OrLandau( lat , &acc , GFINFO.max_iters ,
		    GFINFO.accuracy , Latt.gf_alpha ) ;
#else
  // the memory cheap one wasn't much of a saving so we just use the fast
  iters = Landau( lat , GFINFO.accuracy , GFINFO.max_iters , infile ) ;
#endif
  print_time( ) ;
  return iters ;
}

// print out the details
static void
print_fixing_info( const struct gf_info GFINFO )
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
  case DEFAULT_NOFIX :
    fprintf( stdout , "[GF] No gauge fixing\n" ) ;
    break ;
  case GLU_AXIALT_FIX :
    fprintf( stdout , "[GF] Fixing to lattice temporal gauge.\n"
	     "[GF] WARNING! not a complete gauge fixing prescription.\n" ) ; 
    break ;
  case GLU_COULOMB_FIX :
    fprintf( stdout , "[GF] Coulomb gauge accuracy of %g \n" ,
	     GFINFO.accuracy ) ; 
    break ;
  case GLU_COULOMB_RESIDUAL_FIX :
    fprintf( stdout , "[GF] Coulomb gauge accuracy %g & fixing residual\n" ,
	     GFINFO.accuracy ) ; 
    break ;
  case GLU_LANDAU_FIX :
    fprintf( stdout , "[GF] Landau gauge, accuracy of %g \n" ,
	     GFINFO.accuracy ) ; 
    break ;
  case GLU_MAG_FIX :
    fprintf( stdout , "[GF] Fixing to lattice Maximal Axial Gauge.\n" ) ;
    break ;
  }
  fprintf( stdout , "[GF] Tuning parameter alpha :: %f \n" , Latt.gf_alpha ) ;
  fprintf( stdout , "[GF] Performing AT MOST %zu iterations "
	   "before randomly restarting ... \n" , GFINFO.max_iters ) ;
  fprintf( stdout , "[GF] Allowing for %d restarts before complaint ... \n" , 
	   GF_GLU_FAILURES ) ; 
  // derivative routines available
  fprintf( stdout , "[GF] " ) ;
  #ifdef deriv_lin
  fprintf( stdout , "Using the Hermitian projection definition "
	   "of the gluon fields.\n" ) ; 
  #elif defined deriv_full
  fprintf( stdout , "Using the exact Log definition of the gluon fields.\n" ) ; 
  #endif
  // exponentiation approximation routines
  fprintf( stdout , "[GF] " ) ;
  #if defined exp_approx
  fprintf( stdout , "Approximate O(a) exponential expansion, "
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
  return ;
}

// some simple information for the weird GF routines
static void
output_information( struct site *lat )
{
  fprintf( stdout , "[GF] Tlink %1.15f || Slink %1.15f ||"
	   " Link %1.15f || Plaq %1.15f \n" , 
	   t_links(lat) , s_links(lat) , links(lat) , av_plaquette(lat) ) ;
  return ;
}

// have a wrapper for the above functions here ...
int
GF_wrap( const char *infile , 
	 struct site *lat , 
	 const struct gf_info GFINFO , 
	 const struct head_data HEAD_DATA )
{
  int flag = GLU_SUCCESS ;
  
  print_fixing_info( GFINFO ) ;

  // setup for Landau - Coulomb
  switch( GFINFO.type ) {
  case DEFAULT_NOFIX :
    return GLU_SUCCESS ;
  case GLU_AXIALT_FIX :
    flag = axial_gauge( lat , ND-1 ) ;
    output_information( lat ) ;
    return flag ;
  case GLU_COULOMB_FIX :
    return GF_wrap_coulomb( lat , GFINFO ) ;
  case GLU_COULOMB_RESIDUAL_FIX :
    GF_wrap_coulomb( lat , GFINFO ) ;
    fprintf( stdout , "\n[GF] Residual gauge fixing step \n" ) ;
    residual_fix( lat ) ;
    output_information( lat ) ;
    return GLU_SUCCESS ;
  case GLU_LANDAU_FIX :
    return GF_wrap_landau( lat , infile , GFINFO ) ;
  case GLU_MAG_FIX :
    flag = mag( lat ) ;
    output_information( lat ) ;
    return flag ;
  }
  
  return flag ;
}
