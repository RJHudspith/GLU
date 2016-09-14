/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (GLUlib_wrap.c) is part of GLU.

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
   @file GLUlib_wrap.c
   @brief wrappers for the general functionality of the code
 */
#include "Mainfile.h"

#include "CUT_wrap.h"       // wrap the cutting operations
#include "GLU_memcheck.h"   // basic memory checking
#include "GF_wrap.h"        // wrap the gauge fixing
#include "givens.h"         // allocating the su2 submatrices
#include "par_rng.h"        // initialise_par_rng()
#include "init.h"           // init_latt()
#include "KPHB.h"           // pseudo-heatbath updates
#include "OBS_wrap.h"       // standard observable calculations (default)
#include "random_config.h"  // random transform
#include "read_headers.h"   // read the header data
#include "read_config.h"    // read the binary configuration and check
#include "Scidac.h"         // for the SCIDAC checksums
#include "SM_wrap.h"        // wrap the smearing functions
#include "SUNCxU1_config.h" // compute SU(NC)xU(1) link matrices
#include "taylor_logs.h"    // Taylor series brute force logarithm
#include "writers.h"        // write out a configuration

// allow for some very necessary precomputations ...
#if NC > 3 
  #include "taylor_logs.h"  // precomputes Taylor series coefficients
#endif

// initialisations
void
attach_GLU( void )
{
  // NC-generic precomputations
#if NC > 3
  // factorial computation for the exponentiation routines
  #if !( defined HAVE_LAPACKE_H || defined HAVE_GSL )
  init_factorial( ) ;
  #endif
#endif
  return ;
}

// read a file, has to be out of order because it is called by the others
struct site*
read_file( struct head_data *HEAD_DATA , 
	   const char *config_in )
{
  // some additional information
#ifdef CONDOR_MODE
  fprintf( stdout , "[PAR] CONDOR_MODE selected .... \n" ) ;
#else
  fprintf( stdout , "[PAR] NOT_CONDOR_MODE selected .... "
	   "(same-architecture caching used) \n" ) ;
#endif

  /// here we include the usual stuff look at header for global stuff
  // open our configuration
  FILE *infile = fopen( config_in , "r" ) ;
  /*
  if( infile == NULL ) {
    fprintf( stderr , "[IO] error opening file :: %s\n" , config_in ) ;
    return NULL ;
  }
  */
 
  // initialise the configuration number to zero
  struct head_data tmp ;
  if( read_header( infile , &tmp , GLU_TRUE ) == GLU_FAILURE ) {
    fprintf( stderr , "[IO] Header reading failure \n" ) ;
    fclose( infile ) ;
    return NULL ;
  } 

  // initialise geometry so that we can use LVOLUME and stuff
  init_latt( ) ;

  // check for having enough memory for the gauge field
  if( have_memory_gauge( ) == GLU_FAILURE ) {
    fclose( infile ) ;
    return NULL ;
  }

  // malloc our gauge field and initialise our lattice geometry
  struct site *lat = NULL ;
  if( GLU_malloc( (void**)&lat , 16 , LVOLUME * sizeof( struct site ) ) != 0 ) {
    return NULL ;
  }
  init_navig( lat ) ;

#ifdef SINGLE_PREC
  fprintf( stdout , "[PREC] Single-precision storage for the gauge fields\n" ) ;
#endif

  const int check = get_config_SUNC( infile , lat , tmp ) ;
  // read in the configuration ...  
  if( check == GLU_FAILURE ) {
    fprintf( stderr , "[IO] File read error ... Leaving \n" ) ;
    fclose( infile ) ;
    free( lat ) ;
    return NULL ;
  }

  // look at scidac header again to get the checksums
  // this is taken from the bottom of the file
  if( Latt.head == SCIDAC_HEADER || Latt.head == ILDG_SCIDAC_HEADER ||
      Latt.head == ILDG_BQCD_HEADER ) {
    get_header_data_SCIDAC( infile , &tmp ) ;
  }

  // have a look at some available checks
  if( checks( lat , check , tmp ) == GLU_FAILURE ) { 
    fclose( infile ) ;
    free( lat ) ;
    return NULL ; 
  }

  // set the header info
  *HEAD_DATA = tmp ;

  // free it if it was opened
  if( infile != NULL ) {
    fclose( infile ) ;
  }

  // and finally set the header data into a constant struct
  return lat ;
}

// write out lat
void
write_configuration( struct site *lat , 
		     const char *outfile , 
		     const GLU_output storage , 
		     const char *output_details )
{
  //attempt to open a file
  FILE *out_config = fopen( outfile , "w" ) ;
  write_lat( lat , out_config , storage , output_details ) ;
  fprintf( stdout , "[IO] Configuration written to %s \n", outfile ) ; 
  fclose( out_config ) ;
  return ;
}

// performs heat-bath updates
int
heatbath( const char *infile ,
	  const struct hb_info HBINFO ,
	  const header_mode head , 
	  const char *outfile , 
	  const GLU_output storage , 
	  const char *output_details )
{
  struct head_data HEAD_DATA ;
  struct site *lat = read_file( &HEAD_DATA , infile ) ;

  // should print out a warning
  if( lat == NULL ) return GLU_FAILURE ;

  // heatbath updates
  hb_update( lat , HBINFO , infile , storage , 
	     output_details , (head == UNIT_GAUGE || head == RANDOM_CONFIG) ) ;
  
  // free the gauge fields
  free( lat ) ;

  return GLU_SUCCESS ;
}

// checks unitarity and can write out a configuration
int
read_and_check( const char *infile ,
		const GLU_bool rtrans , 
		const char *outfile , 
		const GLU_output storage , 
		const char *output_details )
{
  struct head_data HEAD_DATA ;
  struct site *lat = read_file( &HEAD_DATA , infile ) ;

  // should print out a warning
  if( lat == NULL ) return GLU_FAILURE ;

  // if we want a random transform then here is where we do it
  if( rtrans == GLU_TRUE ) { random_gtrans( lat ) ; }

  gauge( lat ) ;

  if( Latt.argc == WRITE ) {
    write_configuration( lat , outfile , storage , output_details ) ;
  }

  free( lat ) ;

  return GLU_SUCCESS ;
}

// for the cutting
int
read_and_cut( const char *infile , 
	      const struct cut_info CUTINFO , 
	      const struct sm_info SMINFO )
{
  struct head_data HEAD_DATA ;
  struct site *lat = read_file( &HEAD_DATA , infile ) ;

  // should print out a warning
  if( lat == NULL ) return GLU_FAILURE ;

  cuts_wrap_struct( lat , CUTINFO , SMINFO ) ;

  free( lat ) ;
  return GLU_SUCCESS ;
}

// gauge fixing
int
read_and_fix( const char *infile , 
	      const GLU_bool rtrans , 
	      const struct gf_info GFINFO , 
	      const struct sm_info SMINFO ,
	      const char *outfile , 
	      const GLU_output storage , 
	      const char *output_details )
{
  struct head_data HEAD_DATA ;
  struct site *lat = read_file( &HEAD_DATA , infile ) ;

  if( lat == NULL ) return GLU_FAILURE ;

  // if we want a random transform then here is where we do it
  if( rtrans == GLU_TRUE ) { random_gtrans( lat ) ; }

  GF_wrap( infile , lat , GFINFO , SMINFO , HEAD_DATA ) ;

  if( Latt.argc == WRITE ) {
    write_configuration( lat , outfile , storage , output_details ) ;
  }

  free( lat ) ;
  return GLU_SUCCESS ;
}

// smearing 
int
read_and_smear( const char *infile , 
		const GLU_bool rtrans , 
		const struct sm_info SMINFO ,
		const char *outfile , 
		const GLU_output storage , 
		const char *output_details )
{
  struct head_data HEAD_DATA ;
  struct site *lat = read_file( &HEAD_DATA , infile ) ;

  if( lat == NULL ) return GLU_FAILURE ;

  // if we want a random transform then here is where we do it
  if( rtrans == GLU_TRUE ) { random_gtrans( lat ) ; }

  SM_wrap_struct( lat , SMINFO ) ;

  if( Latt.argc == WRITE ) {
    write_configuration( lat , outfile , storage , output_details ) ;
  }

  free( lat ) ;
  return GLU_SUCCESS ;
}

// U1 the links
int
read_and_U1( const char *infile , 
	     const GLU_bool rtrans ,
	     const struct u1_info U1INFO ,
	     const char *outfile , 
	     const GLU_output storage , 
	     const char *output_details ) 
{
  struct head_data HEAD_DATA ;
  struct site *lat = read_file( &HEAD_DATA , infile ) ;

  if( lat == NULL ) return GLU_FAILURE ;

  // if we want a random transform then here is where we do it
  if( rtrans == GLU_TRUE ) { random_gtrans( lat ) ; }

  suNC_cross_u1( lat , U1INFO ) ;

  if( Latt.argc == WRITE ) {
    write_configuration( lat , outfile , storage , output_details ) ;
  }

  free( lat ) ;
  return GLU_SUCCESS ;
}

// free whatever memory we have allocated
void
unstick_GLU( void )
{
  // free the rng if it has been set
  free_par_rng( ) ; 
  free_latt( ) ;
// if we have to precompute the factorial we free the memory here
#if ( NC > 3 )
  #if !( defined HAVE_LAPACKE_H || defined HAVE_GSL )
  free_factorial( ) ;
  #endif
  free_taylors( ) ;
#endif
  return ;
}
