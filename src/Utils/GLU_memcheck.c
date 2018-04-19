/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (GLU_memcheck.c) is part of GLU.

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
   @file GLU_memcheck.c
   @brief memory checking functions, but at least it is nice to know how much memory our functions will use
 */
#include "Mainfile.h"

// defined in the autoconf procedure
#ifdef HAVE_UNISTD_H
  #include <unistd.h>
  #define DANGEROUS // if we are less than the page size we plough on

  #ifdef GLU_BGQ // BGQ says it has twice as much RAM as it (usually for edinburgh) has ...
  static const double GB = 0.5 * 9.3132257461547852e-10 ; // the mysterious factor is 1.0/(1024*1024*1024)
  #else
  static const double GB = 9.313225746147852e-10 ;
  #endif
static double page_size , MemTotal , MemFree , gauge_fields , gtrans_mats , 
              lat_element , sublat_element ;
static int HYP = -1 , WF = -1 , LGF = -1 , IO = -1 , GFIELDS = -1 , CGF = -1 ;

// call to the OS to request what memory size we have
static void
check_mem( void ) 
{
  const double GLU_size = sizeof( GLU_complex ) ;
  page_size = (double)sysconf( _SC_PAGESIZE ) ; 
  MemTotal = sysconf(_SC_PHYS_PAGES) * page_size ; 
  MemFree = (double)sysconf( _SC_AVPHYS_PAGES ) * page_size ; 
  gauge_fields = GLU_size * (double)LVOLUME * ( NCNC * ND ) * GB + 
    (double)sizeof( int ) * (double)LVOLUME * ( 2 * ND ) * GB ; // integer neighbours
  gtrans_mats = (double)LVOLUME * GLU_size * ( NCNC ) * GB ;
  // lattice sites in GB
  lat_element = GLU_size * (double)LVOLUME * GB ;
  // subcube in GB
  sublat_element = GLU_size * (double)LCU * GB ;
  return ;
}
// and print out the standard message
static void
print_info( const char *info ,
	    const double free_memory ,
	    const double expensive_more ,
	    const double expensive_less )
{
  fprintf( stdout , "[%s] Page size is :: %f, Total malloc "
	   "[More] %f  [Less] %f \n" , 
	   info , ( double )free_memory , expensive_more , expensive_less ) ; 
  return ;
}
#endif

// checks which Landau gauge fixing routine is called
short int 
have_memory_Cgf( void )
{
  int GF = GLU_SUCCESS ;
#ifdef HAVE_UNISTD_H
  if( CGF == -1 ) {
    check_mem( ) ;
    #ifdef DANGEROUS
    const double free_memory = MemTotal * GB ; 
    #else
    const double free_memory = MemFree * GB ;
    #endif
    const double MEM = gauge_fields + ( 3 * NCNC + HERMSIZE ) * sublat_element ;
    if( MEM > free_memory ) {
      fprintf( stdout , "[GF] Avail :: %f    Required :: %f \n" , 
	       free_memory , MEM ) ;
      fprintf( stdout , "[GF] Not enough memory for Coulomb"
	       "gauge fixing ... Leaving\n" ) ;
      GF = GLU_FAILURE ;
    }
  }
#endif
  return GF ;
}

// check we have enough memory
short int 
have_memory_gauge( void )
{
  int gfields = GLU_SUCCESS ;
#ifdef HAVE_UNISTD_H
  if( GFIELDS == -1 ) {
    check_mem( ) ;
    if( gauge_fields > MemTotal * GB ) {
      fprintf( stdout , "[GAUGE] Cannot allocate gauge fields "
	       "[AVAIL] %f [REQUESTED] %f \n" ,
	       MemTotal * GB , gauge_fields ) ;
      gfields = GLU_FAILURE ;
    }
  }
#endif
  return gfields ;
}

// checks which hypercubic blocking code to use ...
short int 
have_memory_hyp( const struct sm_info SMINFO )
{
  int fourDdef = FAST ;
#ifdef HAVE_UNISTD_H 
  if( HYP == -1 ) {
    check_mem( ) ;
    
    //calculate what our most expensive HYP 4D is :: SOMETHING close to this
    const double expensive_less =  (double)NCNC * ND * sizeof( GLU_complex ) * 
      ( LVOLUME * 7./3. + 7. * LCU ) * GB ;
    const double expensive_more =  (double)NCNC * ND * sizeof( GLU_complex ) * 
      ( LVOLUME * 7. + 3. * LCU ) * GB ;
    
    #ifdef DANGEROUS 
    const double free_memory = MemTotal * GB ;
    #else
    const double free_memory = MemFree * GB ;
    #endif	  

    print_info( "SMEAR" , free_memory , expensive_more , expensive_less ) ;
	 
    if( SMINFO.dir == ALL_DIRECTIONS ) {
      if( expensive_more > free_memory ) {
	if( expensive_less > free_memory ) {
	  fprintf( stdout , "[SMEAR] Defined slow smear \n" ) ; 
	  fprintf( stdout , "[SMEAR] Running the REALLY SLOW "
		   "Hypercubically-blocked smearing routine ... \n" ) ; 
	  fourDdef = SLOW ;
	}
        #if NC < 4 // this is not available for NC > 3 because of gratuitous shortening
	else  {
	  fprintf( stdout , "[SMEAR] Running our SLOWER Hypercubically-blocked "
		   "smearing routine ... \n" ) ; 
	  fourDdef = MODERATE ;
	}
	#endif
	// this else if for the expensive_more doodad
      } else {
	fprintf( stdout , "[SMEAR] Running our FASTEST Hypercubically-blocked "
		 "smearing routine ... \n\n" ) ; 
      }
    }
    HYP = fourDdef ;
  } else { fourDdef = HYP ; }
#endif
  return fourDdef ;
}

// checks which Landau gauge fixing routine is called
short int 
have_memory_Lgf( void )
{
  int GF = FAST ;
#ifdef HAVE_UNISTD_H
  if( LGF == -1 ) {
    check_mem( ) ;
    #ifdef DANGEROUS
    const double free_memory = MemTotal * GB ; 
    #else
    const double free_memory = MemFree * GB ;
    #endif

    #if NC == 3
    const int TRUEHERMSIZE = HERMSIZE - 1 ;  
    #else
    const int TRUEHERMSIZE = HERMSIZE ;
    #endif

    //calculate what our most expensive Gauge fixer is
    #ifdef LUXURY_GAUGE
    const double expensive_less = 3*gauge_fields + gtrans_mats + lat_element * 2.0 ;
    const double expensive_more = 3*gauge_fields + gtrans_mats + TRUEHERMSIZE * lat_element * 2.0 ;
    #else
    const double expensive_less = gauge_fields + gtrans_mats + lat_element * 2.0 ;
    const double expensive_more = gauge_fields + gtrans_mats + TRUEHERMSIZE * lat_element * 2.0 ;
    #endif

    print_info( "GF" , free_memory , expensive_more , expensive_less ) ;
  
    if ( expensive_more < free_memory ) {
      fprintf( stdout , "[GF] Running our fastest Landau Gauge Fixing "
	       "routine .. \n\n" ) ; 
    } else {
      fprintf( stdout , "[GF] Sorry, not enough memory to run Landau Gauge "
	       "Fixing ... \n\n" ) ;
      GF = GLU_FAILURE ;
    }
    LGF = GF ;
  } else { GF = LGF ; } 
#endif
  return GF ;
} 

// checks which Landau gauge fixing routine is called
short int 
have_memory_readers_writers( const GLU_output config_type )
{
  int CONFIG = FAST ;
#ifdef HAVE_UNISTD_H
  if( IO == -1 ) {
    check_mem( ) ;
    
    double expensive_more , expensive_less ;
    
    #ifdef DANGEROUS
    const double free_memory = MemTotal * GB ;
    #else
    const double free_memory = MemFree * GB ;
    #endif

    // check to make sure we can fit gauge field and the malloc
    IO = FAST ;
    switch( config_type ) 
      {
      case OUTPUT_SMALL :
	expensive_more = gauge_fields + (double)sizeof( double ) * ( (double)LVOLUME * ( ND * ( NCNC - 1 ) ) ) * GB ;
	expensive_less = gauge_fields + (double)sizeof( double ) * ( (double)( ND * ( NCNC - 1 ) ) ) * GB ;
	print_info( "IO" , free_memory , expensive_more , expensive_less ) ;
	if( expensive_more < free_memory ) {
	  return FAST ;
	} else { 
	  IO = SLOW ;
	  return SLOW ;
	}
      case OUTPUT_GAUGE :
	expensive_more = gauge_fields + (double)sizeof( double complex ) * ( (double)LVOLUME * ( ND * NC * ( NC - 1 ) ) ) * GB ;
	expensive_less = gauge_fields + (double)sizeof( double complex ) * ( (double)( ND * NC * ( NC - 1 ) ) ) * GB ;
	print_info( "IO" , free_memory , expensive_more , expensive_less ) ;
	if( expensive_more < free_memory ) {
	  return FAST ;
	} else { 
	  IO = SLOW ;
	  return SLOW ;
	}
      case OUTPUT_MILC :
      case OUTPUT_SCIDAC :
      case OUTPUT_ILDG :
      case OUTPUT_NCxNC :
	expensive_more = 2.0 * gauge_fields ; 
	expensive_less = gauge_fields + (double) sizeof( double complex ) * ( (double)( ND * NCNC ) ) * GB ;
	print_info( "IO" , free_memory , expensive_more , expensive_less ) ;
	if( expensive_more < free_memory ) {
	  return FAST ;
	} else { 
	  IO = SLOW ;
	  return SLOW ;
	}
      case OUTPUT_HIREP :
	return FAST ;
      default :
	return SLOW ;
      } 
  } else { CONFIG = IO ; }
#endif
  return CONFIG ;
} 

// checks which wilson flow routine is called
short int 
have_memory_wf( const struct sm_info SMINFO ) 
{
  int fourDdef = FAST ; 
#ifdef HAVE_UNISTD_H
  if( WF == -1 ) {
    check_mem( ) ;

    #if NC == 3 
    const double TRUEHERM = HERMSIZE - 1 ;
    #else
    const double TRUEHERM = HERMSIZE ;
    #endif

    double expensive_less = gauge_fields + ( TRUEHERM ) * lat_element * ND + ( 3.0 * sublat_element * NCNC * ND ) ;
    double expensive_more = 2.0 * gauge_fields + TRUEHERM * lat_element * ND ;
    
    if( SMINFO.type == SM_ADAPTWFLOW_LOG || SMINFO.type == SM_ADAPTWFLOW_STOUT ) {
      expensive_less += gauge_fields ;
      expensive_more += gauge_fields ;
    }

    #ifdef DANGEROUS
    const double free_memory = MemTotal * GB ;
    #else
    const double free_memory = Memfree * GB ;
    #endif

    print_info( "WFLOW" , free_memory , expensive_more , expensive_less ) ;
    if( SMINFO.dir == ALL_DIRECTIONS )  {
      if( expensive_more > free_memory ) {
	if( expensive_less > free_memory ) {
	  fprintf( stdout , "[WFLOW] Cannot use either Wilson Flow "
		   "routine .. \n\n" ) ; 
	  fourDdef = GLU_FAILURE ;
	} else {
	  fprintf( stdout , "[WFLOW] Running our slower Wilson Flow "
		   "routine .. \n\n" ) ; 
	  fourDdef = MODERATE ;
	}
      } else {
	fprintf( stdout , "[WFLOW] Running our fastest Wilson Flow "
		 "routine .. \n\n" ) ; 
      }
    }
    WF = fourDdef ;
  } else { fourDdef = WF ; } 
#endif
  return fourDdef ;
}

// clean this up for local scoping
#ifdef DANGEROUS
  #undef DANGEROUS
#endif
