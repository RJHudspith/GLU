/*
    Copyright 2013 Renwick James Hudspith

    This file (input_reader.c) is part of GLU.

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
   @file input_reader.c
   @brief reader for the input file

   This file is too long, I should split it up
*/

#include "Mainfile.h"

#include <errno.h>
#include "GLU_rng.h"  // we initialise the RNG from the input file

// we might want to change this at some point
#define GLU_STR_LENGTH 64

// maximum number of tokens
#define INPUTS_LENGTH 36

// tokenize the input file
struct inputs {
  char TOKEN[ GLU_STR_LENGTH ] ;
  char VALUE[ GLU_STR_LENGTH ] ;
} ;

// is this what valgrind dislikes?
static struct inputs *INPUT = NULL ;

// counter for the number of tags
static int NTAGS = 0 ;

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// allocate the input file struct
static void 
pack_inputs( FILE *setup )
{
  INPUT = ( struct inputs* )calloc( INPUTS_LENGTH , sizeof( struct inputs ) ) ;
  // and put into the buffer
  while( NTAGS++ , fscanf( setup , "%s = %s" , INPUT[ NTAGS ].TOKEN , 
			   INPUT[ NTAGS ].VALUE )  != EOF ) { }
  return ;
}

// frees up the struct
static void
unpack_inputs( void )
{
  free( INPUT ) ;
  return ;
}

// prints out the problematic tag
static int
tag_failure( const char *tag )
{
  printf( "[IO] Failure looking for tag %s in input file ... Leaving\n" , tag ) ;
  return GLU_FAILURE ;
}

// look for a tag and return the index if found, else return -1
static int
tag_search( const char *tag ) 
{
  int i ;
  for( i = 0 ; i < NTAGS ; i++ ) {
    if( are_equal( INPUT[i].TOKEN , tag ) ) return i ;
  }
  return GLU_FAILURE ;
}

// I use this pattern quite a bit
static int
setint( size_t *val ,
	const char *tag )
{
  char *endptr ;
  errno = 0 ;
  const int idx = tag_search( tag ) ;
  if( idx == GLU_FAILURE ) { return tag_failure( tag ) ; }
  *val = (size_t)strtol( INPUT[idx].VALUE , &endptr , 10 ) ;
  if( endptr == INPUT[idx].VALUE || errno == ERANGE ) {
    return GLU_FAILURE ;
  }
  return GLU_SUCCESS ;
}

// above but for doubles
static int
setdbl( double *val ,
	const char *tag )
{
  char *endptr ;
  errno = 0 ;
  const int idx = tag_search( tag ) ;
  if( idx == GLU_FAILURE ) { return tag_failure( tag ) ; }
  *val = strtod( INPUT[idx].VALUE , &endptr ) ;
  if( endptr == INPUT[idx].VALUE || errno == ERANGE ) {
    return GLU_FAILURE ;
  }
  return GLU_SUCCESS ;
}

// quickly get the configuration number from the input file
// must be greater than or equal to zero
static size_t
confno( void )
{
  size_t conf = 0 ;
  if( setint( &conf , "CONFNO" ) == GLU_FAILURE ) {
    printf( "[IO] I do not understand CONFNO %zu \n" , conf ) ;
    printf( "[IO] CONFNO should be greater than 0\n" ) ;
    return GLU_FAILURE ;
  }
  return conf ;
}

// read from the input file our dimensions
static int
read_random_lattice_info( void )
{
  // lattice dimensions set from input file
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    char tmp[ GLU_STR_LENGTH ] ;
    sprintf( tmp , "DIM_%d" , mu ) ;
    Latt.dims[mu] = 8 ; // default to an 8^{ND} lattice
    if( setint( &Latt.dims[ mu ] , tmp ) == GLU_FAILURE ) {
      return GLU_FAILURE ;
    }
  }
  if( setint( &Latt.flow , "CONFNO" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  return GLU_SUCCESS ;
}

//////////////////// input file functions /////////////////////

// get the information for the header
static int
config_information( char *details )
{
  const int details_idx = tag_search( "CONFIG_INFO" ) ;
  if( are_equal( INPUT[details_idx].VALUE , "YES" ) ) { 
    return tag_failure( "CONFIG_INFO" ) ; 
  }
  sprintf( details , "%s" , INPUT[details_idx].VALUE ) ;
  return GLU_SUCCESS ;
}

// detect what mode we are using
static int
get_mode( GLU_mode *mode )
{
  // look for which mode we are running in
  {
    const int mode_idx = tag_search( "MODE" ) ;
    if( mode_idx == GLU_FAILURE ) { return tag_failure( "MODE" ) ; }
    // what mode do we want?
    if( are_equal( INPUT[mode_idx].VALUE , "GAUGE_FIXING" ) ) {
      *mode = MODE_GF ;
    } else if( are_equal( INPUT[mode_idx].VALUE , "CUTTING" ) ) {
      *mode = MODE_CUTS ;
    } else if( are_equal( INPUT[mode_idx].VALUE , "SMEARING" ) ) {
      *mode = MODE_SMEARING ;
    } else if( are_equal( INPUT[mode_idx].VALUE , "SUNCxU1" ) ) {
      *mode = MODE_CROSS_U1 ;
    } else {
      *mode = MODE_REWRITE ;
    }
  }
  // look at what seed type we are going to use
  {
    const int seed_idx = tag_search( "SEED" ) ;
    if( seed_idx == GLU_FAILURE ) { return tag_failure( "SEED" ) ; }
    sscanf( INPUT[seed_idx].VALUE , "%u" , &Latt.Seed[0] ) ;
    if( Latt.Seed[0] == 0 ) {
      printf( "\n[RNG] Generating RNG seed from urandom \n" ) ; 
      if( initialise_seed( ) == GLU_FAILURE ) return GLU_FAILURE ;
    }
    #ifdef KISS_RNG
    printf( "[RNG] KISS Seed %u \n\n" , Latt.Seed[0] ) ;
    #elif defined MWC_1038_RNG
    printf( "[RNG] MWC_1038 Seed %u \n\n" , Latt.Seed[0] ) ;
    #elif defined MWC_4096_RNG
    printf( "[RNG] MWC_4096 Seed %u \n\n" , Latt.Seed[0] ) ;
    #elif defined GSL_RNG
    printf( "[RNG] GSL (MT) Seed %u \n\n" , Latt.Seed[0] ) ;
    #elif defined XOR1024_RNG
    printf( "[RNG] XOR1024 Seed %u \n\n" , Latt.Seed[0] ) ;
    #else
    printf( "[RNG] well_19937a Seed %u \n\n" , Latt.Seed[0] ) ;
    #endif
  }
  return GLU_SUCCESS ;
}

// the code that writes out all of the config details
static GLU_output
out_details( const GLU_mode mode )
{  
  // get the storage type
  {
    const int storage_idx = tag_search( "STORAGE" ) ;
    if( storage_idx == GLU_FAILURE ) { return tag_failure( "STORAGE" ) ; }
    if( are_equal( INPUT[storage_idx].VALUE , "NERSC_SMALL" ) ) {
      // small is not available for larger NC default to NCxNC
      #if NC > 3
      return OUTPUT_NCxNC ;
      #else
      return ( mode != MODE_CROSS_U1 ) ? OUTPUT_SMALL : OUTPUT_NCxNC ;
      #endif
    } else if( are_equal( INPUT[storage_idx].VALUE , "NERSC_GAUGE" ) ) {
      return ( mode != MODE_CROSS_U1 ) ? OUTPUT_GAUGE : OUTPUT_NCxNC ;
    } else if( are_equal( INPUT[storage_idx].VALUE , "NERSC_NCxNC" ) ) {
      return OUTPUT_NCxNC ;
    } else if( are_equal( INPUT[storage_idx].VALUE , "HIREP" ) ) {
      return OUTPUT_HIREP ;
    } else if( are_equal( INPUT[storage_idx].VALUE , "MILC" ) ) {
      return OUTPUT_MILC ;
    } else if( are_equal( INPUT[storage_idx].VALUE , "SCIDAC" ) ) {
      return OUTPUT_SCIDAC ;
    } else if( are_equal( INPUT[storage_idx].VALUE , "ILDG" ) ) {
      return OUTPUT_ILDG ;
    } else {
      // default to NERSC_NCxNC
      return OUTPUT_NCxNC ;
    }
  }
  // we will not reach this point
  return GLU_FAILURE ;
}

// get the header type
static header_mode
header_type( void )
{
  // headers we support
  {
    const int header_idx = tag_search( "HEADER" ) ;
    if( header_idx == GLU_FAILURE ) { return tag_failure( "HEADER" ) ; }
    if( are_equal( INPUT[header_idx].VALUE , "NERSC" ) ) {
      printf( "[IO] Attempting to read a NERSC file \n" ) ;
      return NERSC_HEADER ;
    } else if( are_equal( INPUT[header_idx].VALUE , "HIREP" ) ) {
      if( ( Latt.flow = confno( ) ) == GLU_FAILURE ) return UNSUPPORTED ;
      printf( "[IO] Attempting to read a HIREP file \n" ) ;
      printf( "[IO] Using sequence number from input file :: %zu \n" ,
	      Latt.flow ) ;
      return HIREP_HEADER ;
    } else if( are_equal( INPUT[header_idx].VALUE , "MILC" ) ) {
      if( ( Latt.flow = confno( ) ) == GLU_FAILURE ) return UNSUPPORTED ;
      printf( "[IO] Attempting to read a MILC file \n" ) ;
      printf( "[IO] Using sequence number from input file :: %zu \n" ,
	      Latt.flow ) ;
      return MILC_HEADER ;
    } else if( are_equal( INPUT[header_idx].VALUE , "SCIDAC" ) ) {
      if( ( Latt.flow = confno( ) ) == GLU_FAILURE ) return UNSUPPORTED ;
      printf( "[IO] Attempting to read a SCIDAC file \n" ) ;
      printf( "[IO] Using sequence number from input file :: %zu \n" ,
	      Latt.flow ) ;
      return SCIDAC_HEADER ;
    } else if( are_equal( INPUT[header_idx].VALUE , "LIME" ) ) {
      if( ( Latt.flow = confno( ) ) == GLU_FAILURE ) return UNSUPPORTED ;
      printf( "[IO] Attempting to read an LIME file \n" ) ;
      printf( "[IO] Using sequence number from input file :: %zu \n" ,
	      Latt.flow ) ;
      printf( "[IO] WARNING!! NOT CHECKING ANY CHECKSUMS!! \n" ) ;
      return LIME_HEADER ;
    } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_SCIDAC" ) ) {
      if( ( Latt.flow = confno( ) ) == GLU_FAILURE ) return UNSUPPORTED ;
      printf( "[IO] Attempting to read an ILDG (Scidac) file \n" ) ;
      printf( "[IO] Using sequence number from input file :: %zu \n" ,
	      Latt.flow ) ;
      return ILDG_SCIDAC_HEADER ;
    } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_BQCD" ) ) {
      if( ( Latt.flow = confno( ) ) == GLU_FAILURE ) return UNSUPPORTED ;
      printf( "[IO] Attempting to read an ILDG (BQCD) file \n" ) ;
      printf( "[IO] Using sequence number from input file :: %zu \n" ,
	      Latt.flow ) ;
      return ILDG_BQCD_HEADER ;
    } else if( are_equal( INPUT[header_idx].VALUE , "RANDOM" ) ) {
      printf( "[IO] Attempting to generate an SU(%d) RANDOM config \n" , NC ) ;
      if( read_random_lattice_info( ) == GLU_FAILURE ) return UNSUPPORTED ;
      return RANDOM_CONFIG ;
    } else if( are_equal( INPUT[header_idx].VALUE , "UNIT" ) ) {
      printf( "[IO] Attempting to generate an %dx%d UNIT config \n" , NC , NC ) ;
      read_random_lattice_info( ) ;
      return UNIT_GAUGE ;
    } else if( are_equal( INPUT[header_idx].VALUE , "INSTANTON" ) ) {
      printf( "[IO] Attempting to generate a SU(%d) BPST instanton config \n" 
	      , NC ) ;
      if( read_random_lattice_info( ) == GLU_FAILURE ) return UNSUPPORTED ;
      return INSTANTON ;
    }
    printf( "[IO] HEADER %s not recognised ... Leaving \n" , 
	    INPUT[header_idx].VALUE ) ;
    return UNSUPPORTED ;
  }
  return UNSUPPORTED ; 
}

// pack the cut_info struct
static int
read_cuts_struct( struct cut_info *CUTINFO )
{
  // set up the cuttype
  {
    const int cuttype_idx = tag_search( "CUTTYPE" ) ;
    if( cuttype_idx == GLU_FAILURE ) { return tag_failure( "CUTTYPE" ) ; }
    if( are_equal( INPUT[cuttype_idx].VALUE , "EXCEPTIONAL" ) ) {
      CUTINFO -> dir = EXCEPTIONAL ;
    } else if( are_equal( INPUT[cuttype_idx].VALUE , "NON_EXCEPTIONAL" ) ) {
      CUTINFO -> dir = NONEXCEPTIONAL ;
    } else if( are_equal( INPUT[cuttype_idx].VALUE , "FIELDS" ) ) {
      CUTINFO -> dir = FIELDS ;
    } else if( are_equal( INPUT[cuttype_idx].VALUE , "SMEARED_GLUONS" ) ) {
      CUTINFO -> dir = SMEARED_GLUONS ;
    } else if( are_equal( INPUT[cuttype_idx].VALUE , "INSTANTANEOUS_GLUONS" ) ) {
      CUTINFO -> dir = INSTANTANEOUS_GLUONS ;
    } else if( are_equal( INPUT[cuttype_idx].VALUE , "CONFIGSPACE_GLUONS" ) ) {
      CUTINFO -> dir = CONFIGSPACE_GLUONS ;
    } else if( are_equal( INPUT[cuttype_idx].VALUE , "GLUON_PROPS" ) ) {
      CUTINFO -> dir = GLUON_PROPS ;
    } else if( are_equal( INPUT[cuttype_idx].VALUE , "STATIC_POTENTIAL" ) ) {
      CUTINFO -> dir = STATIC_POTENTIAL ;
    } else if( are_equal( INPUT[cuttype_idx].VALUE , "TOPOLOGICAL_SUSCEPTIBILITY" ) ) {
      CUTINFO -> dir = TOPOLOGICAL_SUSCEPTIBILITY ;
    } else {
      printf( "[IO] I do not understand your CUTTYPE %s\n" , 
	      INPUT[cuttype_idx].VALUE ) ;
      printf( "[IO] Defaulting to no cutting \n" ) ;
      CUTINFO -> dir = GLU_FAILURE ;
    }
  }
  // momentum space cut def
  {
    const int momcut_idx = tag_search( "MOM_CUT" ) ;
    if( momcut_idx == GLU_FAILURE ) { return tag_failure( "MOM_CUT" ) ; }
    if ( are_equal( INPUT[momcut_idx].VALUE , "HYPERCUBIC_CUT" ) ) {
      CUTINFO -> type = HYPERCUBIC_CUT ; 
    } else if ( are_equal( INPUT[momcut_idx].VALUE , "SPHERICAL_CUT" ) ) {
      CUTINFO -> type = PSQ_CUT ; 
    } else if ( are_equal( INPUT[momcut_idx].VALUE , "CYLINDER_CUT" ) ) {
      CUTINFO -> type = CYLINDER_CUT ; 
    } else if ( are_equal( INPUT[momcut_idx].VALUE , "CONICAL_CUT" ) ) {
      CUTINFO -> type = CYLINDER_AND_CONICAL_CUT ; 
    } else {
      printf( "[IO] Unrecognised type [%s] \n" , INPUT[momcut_idx].VALUE ) ; 
      printf( "[IO] Defaulting to SPHERICAL_CUT \n" ) ; 
    }
  }
  // field definition
  {
    const int field_idx = tag_search( "FIELD_DEFINITION" ) ;
    if( field_idx == GLU_FAILURE ) { return tag_failure( "FIELD_DEFINITION" ) ; }
    CUTINFO -> definition = LINEAR_DEF ;
    if( are_equal( INPUT[field_idx].VALUE , "LOGARITHMIC" ) ) {
      CUTINFO -> definition = LOG_DEF ;
    } 
  }
  // minmom, maxmom angle and cylinder width
  if( setint( &( CUTINFO -> max_t ) , "MAX_T" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  if( setint( &( CUTINFO -> max_mom ) , "MAXMOM" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  if( setint( &( CUTINFO -> angle ) , "ANGLE" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  // set the cylinder width
  if( setdbl( &( CUTINFO -> cyl_width ) , "CYL_WIDTH" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  // look for where the output is going
  {
    const int output_idx = tag_search( "OUTPUT" ) ;
    if( output_idx == GLU_FAILURE ) { return tag_failure( "OUTPUT" ) ; }
    sprintf( CUTINFO -> where , "%s" , INPUT[output_idx].VALUE ) ; 
  }
  return GLU_SUCCESS ;
}

// read the gauge fixing information
static int
read_gf_struct ( struct gf_info *GFINFO )
{
  // look at what seed type we are going to use
  {
    const int gf_idx = tag_search( "GFTYPE" ) ;
    if( gf_idx == GLU_FAILURE ) { return tag_failure( "GFTYPE" ) ; } 
    if( are_equal( INPUT[gf_idx].VALUE , "LANDAU" ) ) { 
      GFINFO -> type = GLU_LANDAU_FIX ;
    } else if ( are_equal( INPUT[gf_idx].VALUE , "COULOMB" ) ) {
      GFINFO -> type = GLU_COULOMB_FIX ;
    } else {
      printf( "[IO] unknown type [%s] : Defaulting to NO GAUGE FIXING \n" , 
	      INPUT[gf_idx].VALUE ) ; 
      GFINFO -> type = DEFAULT_NOFIX ; 
    }
  }
  // have a look to see what "improvements" we would like
  {
    const int improve_idx = tag_search( "IMPROVEMENTS" ) ;
    if( improve_idx == GLU_FAILURE ) { return tag_failure( "IMPROVEMENTS" ) ; } 
    GFINFO -> improve = NO_IMPROVE ; // default is no "Improvements" 
    if( are_equal( INPUT[improve_idx].VALUE , "MAG" ) ) {
      GFINFO -> improve = MAG_IMPROVE ; 
    } else if ( are_equal( INPUT[improve_idx].VALUE , "SMEAR" ) ) {
      GFINFO -> improve = SMPREC_IMPROVE ; 
    } else if(  are_equal( INPUT[improve_idx].VALUE , "RESIDUAL" ) ) {
      GFINFO -> improve = RESIDUAL_IMPROVE ; 
    }
  }
  // set the accuracy is 10^{-ACCURACY}
  if( setdbl( &( GFINFO -> accuracy ) , "ACCURACY" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  } 
  GFINFO -> accuracy = pow( 10.0 , -GFINFO -> accuracy ) ;
  // set the maximum number of iterations of the routine
  if( setint( &( GFINFO -> max_iters ) , "MAX_ITERS" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  // set the alpha goes in Latt.gf_alpha for some reason
  if( setdbl( &Latt.gf_alpha , "GF_TUNE" ) == GLU_FAILURE ) {
    printf( "Failure \n" ) ;
    return GLU_FAILURE ;
  }
  return GLU_SUCCESS ;
}

// the one for the U(1) data
static int
read_suNC_x_U1( struct u1_info *U1INFO )
{
  // U1 coupling strength
  if( setdbl( &( U1INFO -> alpha ) , "U1_ALPHA" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  // U1 charge
  if( setdbl( &( U1INFO -> charge ) , "U1_CHARGE" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  // U1 measurement type default is just the plaquette
  {
    const int U1_meas_idx = tag_search( "U1_MEAS" ) ;
    if( U1_meas_idx == GLU_FAILURE ) { return tag_failure( "U1_MEAS" ) ; }
    U1INFO -> meas = U1_PLAQUETTE ;
    if( are_equal( INPUT[ U1_meas_idx ].VALUE , "U1_RECTANGLE" ) ) {
      U1INFO -> meas = U1_RECTANGLE ;
    } else if ( are_equal( INPUT[ U1_meas_idx ].VALUE , "U1_TOPOLOGICAL" ) ) {
      U1INFO -> meas = U1_TOPOLOGICAL ;
    }
  }
  return GLU_SUCCESS ;
}

// are we performing a random transform?
static GLU_bool
rtrans( void )
{
  const int rtrans_idx = tag_search( "RANDOM_TRANSFORM" ) ;
  if( are_equal( INPUT[rtrans_idx].VALUE , "YES" ) ) return GLU_TRUE ;
  return GLU_FALSE ;
}

// pack the sm_info struct
static int
smearing_info( struct sm_info *SMINFO )
{
  // find the smeartype index
  {
    const int type_idx = tag_search( "SMEARTYPE" ) ;
    if( type_idx == GLU_FAILURE ) { return tag_failure( "SMEARTYPE" ) ; }
    if( are_equal( INPUT[type_idx].VALUE , "APE" ) ) { 
      SMINFO -> type = SM_APE ;
    } else if( are_equal( INPUT[type_idx].VALUE , "STOUT" ) ) {
      SMINFO -> type = SM_STOUT ;
    } else if( are_equal( INPUT[type_idx].VALUE , "LOG" ) ) {
      SMINFO -> type = SM_LOG ;
    } else if( are_equal( INPUT[type_idx].VALUE , "HYP" ) ) {
      SMINFO -> type = SM_HYP ;
    } else if( are_equal( INPUT[type_idx].VALUE , "HEX" ) ) {
      SMINFO -> type = SM_HEX ;
    } else if( are_equal( INPUT[type_idx].VALUE , "HYL" ) ) {
      SMINFO -> type = SM_HYL ;
    } else if( are_equal( INPUT[type_idx].VALUE , "WFLOW_LOG" ) ) {
      SMINFO -> type = SM_WFLOW_LOG ;
    } else if( are_equal( INPUT[type_idx].VALUE , "WFLOW_STOUT" ) ) {
      SMINFO -> type = SM_WFLOW_STOUT ;
    } else if( are_equal( INPUT[type_idx].VALUE , "ADAPTWFLOW_LOG" ) ) {
      SMINFO -> type = SM_ADAPTWFLOW_LOG ;
    } else if( are_equal( INPUT[type_idx].VALUE , "ADAPTWFLOW_STOUT" ) ) {
      SMINFO -> type = SM_ADAPTWFLOW_STOUT ;
    } else {
      printf( "[IO] Unrecognised Type [%s] Defaulting to No Smearing \n" , 
	      INPUT[type_idx].VALUE ) ;
      SMINFO -> type = SM_NOSMEARING ;
    }
  }
  // look for the number of directions
  {
    const int dir_idx = tag_search( "DIRECTION" ) ;
    if( dir_idx == GLU_FAILURE ) { return tag_failure( "DIRECTION" ) ; }
    if( are_equal( INPUT[dir_idx].VALUE , "SPATIAL" ) ) {
      SMINFO -> dir = SPATIAL_LINKS_ONLY ;
    } else {
      SMINFO -> dir = ALL_DIRECTIONS ;
    }
  }
  // set up the number of smearing iterations
  if( setint( &( SMINFO -> smiters ) , "SMITERS" ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  // poke the smearing alpha's into Latt.smalpha[ND]
  // logically for an ND - dimensional theory there are ND - 1 HYP params.
  {
    int mu ;
    for( mu = 0 ; mu < ND - 1 ; mu++ ) {
      char alpha_str[ 64 ] ;
      sprintf( alpha_str , "ALPHA%d" , mu + 1 ) ;
      if( setdbl( &Latt.sm_alpha[ mu ] , alpha_str ) == GLU_FAILURE ) {
	return GLU_FAILURE ;
      }
    }
  }
  return GLU_SUCCESS ;
}

// fills the INFILE struct with all the useful information
int
get_input_data( struct infile_data *INFILE ,
		const char *file_name )
{
  // open the input file in here and free it at the bottom
  FILE *infile = fopen( file_name , "r" ) ;
  if( infile == NULL ) {
    printf( "[IO] input file cannot be read ... Leaving\n" ) ;
    return GLU_FAILURE ;
  }

  // if we can open the file we push it into a big structure
  pack_inputs( infile ) ;

  int INPUT_FAILS = 0 ; // counter for the number of failures

  // fill in the input data into a struct
  if( get_mode( &( INFILE -> mode ) ) == GLU_FAILURE ) INPUT_FAILS++ ;

  if( read_gf_struct ( &( INFILE -> GFINFO ) ) == GLU_FAILURE ) INPUT_FAILS++ ;

  if( read_cuts_struct( &( INFILE -> CUTINFO ) ) == GLU_FAILURE ) INPUT_FAILS++ ;

  if( read_suNC_x_U1( &( INFILE -> U1INFO ) ) == GLU_FAILURE ) INPUT_FAILS++ ;

  if( smearing_info( &( INFILE -> SMINFO ) ) == GLU_FAILURE ) INPUT_FAILS++ ;
 
  // are we performing a random transform
  INFILE -> rtrans = rtrans( ) ;

  // get the header type
  INFILE -> head = header_type( ) ;
  if( INFILE -> head == UNSUPPORTED ) INPUT_FAILS++ ;
  Latt.head = INFILE -> head ; // set the header type

  // and the storage type
  INFILE -> storage = out_details( INFILE -> mode ) ;

  // put the config info 
  config_information( INFILE -> output_details ) ;

  // close the file and deallocate the buffer
  fclose( infile ) ;
  unpack_inputs( ) ;

  // if we have hit ANY problem we return GLU_FAILURE this causes it to exit
  return ( INPUT_FAILS != 0 ) ? GLU_FAILURE : GLU_SUCCESS ;
}

#undef GLU_STR_LENGTH
#undef INPUTS_LENGTH
