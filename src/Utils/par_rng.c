/**
   @file par_rng.c
   @brief parallel RNG codes live here
   @warning these are parallelised the stupid way

   TODO :: write it in such a way that if we don't have openmp this still works
 */
#include "Mainfile.h"

#include "chklat_stuff.h" // reading a NERSC header
#include "GLU_bswap.h"    // for byte swapping
#include "par_rng.h"
#include "par_KISS.h" // parallel KISS generator
#include "par_MWC_1038.h" // parallel MWC generator
#include "par_MWC_4096.h" // parallel MWC generator
#include "par_XOR_1024.h" // parallel XOR generator
#include "par_WELL_512.h" // parallel well rng

// have we initialised the rng?
static GLU_bool RNG_inited = GLU_FALSE ;

// seed the rng
int
initialise_par_rng( const char *rng_file ) 
{
  if( RNG_inited == GLU_FALSE ) {

    // tell us our rng
#if (defined KISS_RNG)
    fprintf( stdout , "[PAR_RNG] KISS_RNG \n" ) ;
#elif (defined MWC_4096_RNG)
    fprintf( stdout , "[PAR_RNG] MWC_4096\n" ) ;
#elif (defined MWC_1038_RNG)
    fprintf( stdout , "[PAR_RNG] MWC_1038\n" ) ;
#elif (defined XOR_1024_RNG)
    fprintf( stdout , "[PAR_RNG] XOR_1024\n" ) ;
#else
    fprintf( stdout , "[PAR_RNG] WELL_512\n" ) ;
#endif

    if( rng_file == NULL ) {
      // pull from the entropy pool?
      uint32_t *Seeds = malloc( Latt.Nthreads * sizeof( uint32_t ) ) ;

      size_t i ;
      if( Latt.Seed[0] == 0 ) {
	FILE *urandom = fopen( "/dev/urandom" , "r" ) ;
	if( urandom == NULL ) {
	  fprintf( stderr , "[RNG] /dev/urandom not opened!! ... Exiting \n" ) ;
	  return GLU_FAILURE ;
	}
	// read them from urandom
	if( fread( Seeds , sizeof( Latt.Seed ) , 
		   Latt.Nthreads , urandom ) != Latt.Nthreads ) { 
	  fprintf( stderr , "[RNG] Entropy pool Seed not read properly ! "
		   "... Exiting \n" ) ; 
	  return GLU_FAILURE ;
	}
	fclose( urandom ) ;
      } else {
	//
	for( i = 0 ; i < Latt.Nthreads ; i++ ) {
	  Seeds[ i ] = Latt.Seed[0] + i ;
	}
	//
      }

      fprintf( stdout , "[PAR_RNG] Entropy read \n" ) ;
      for( i = 0 ; i < Latt.Nthreads ; i++ ) {
	fprintf( stdout , "[PAR_RNG] Seed_%zu %u \n" , i , Seeds[i] ) ;
      }

      // do the seeding
      #if (defined KISS_RNG)
      GLU_set_par_KISS_table( Seeds ) ;
      #elif (defined MWC_4096_RNG)
      GLU_set_par_MWC_4096_table( Seeds ) ;
      #elif (defined MWC_1038_RNG)
      GLU_set_par_MWC_1038_table( Seeds ) ;
      #elif (defined XOR_1024_RNG)
      GLU_set_par_XOR_1024_table( Seeds ) ;
      #else
      GLU_set_par_WELL_512_table( Seeds ) ;
      #endif

      // warm up the rng
      #pragma omp parallel for private(i)
      for( i = 0 ; i < Latt.Nthreads ; i++ ) {
	size_t j ;
	const uint32_t thread = get_GLU_thread( ) ;
	for( j = 0 ; j < 10000 ; j++ ) {
	  par_rng_dbl( thread ) ;
	}
      }
      fprintf( stdout , "[PAR_RNG] warmed up\n" ) ;

    } else {
      return read_par_rng_state( rng_file ) ;
    }
    RNG_inited = GLU_TRUE ;
  }

  return GLU_SUCCESS ;
}

// free the table
void
free_par_rng( void )
{
  if( RNG_inited == GLU_TRUE ) {
#if (defined KISS_RNG)
    free_par_KISS( ) ;
#elif (defined MWC_4096_RNG)
    free_par_MWC_4096( ) ;
#elif (defined XOR_1024_RNG)
    free_par_XOR_1024( ) ;
#elif (defined MWC_1038_RNG)
    free_par_MWC_1038( ) ;
#else
    free_par_WELL_512( ) ;
#endif
  }
  return ;
}

// generate a random NCxNC matrix ... Use gaussian numbers?
void 
par_generate_NCxNC( GLU_complex U[ NCNC ] , 
		    const uint32_t thread )
{
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    // gaussians
    U[i] = par_polar_box( thread ) ;
  }
  return ;
}

// Gaussian distributed doubles ocassionally called
GLU_real
par_polar( const uint32_t thread )
{
  const GLU_real u = par_rng_dbl( thread ) , v = par_rng_dbl( thread ) ;
  return sqrt( -2. * log( u ) ) * cos( TWOPI * v ) ;
}

// looks more complicated , is faster.
GLU_complex
par_polar_box( const uint32_t thread )
{ 
  register const GLU_real u = (GLU_real)( 2. * par_rng_dbl( thread ) - 1. ) ;
  register const GLU_real v = (GLU_real)( 2. * par_rng_dbl( thread ) - 1. ) ;
  const GLU_real s = u * u + v * v ;
  return s < 1. ? sqrt( -log( s ) / s ) * ( u + I * v ) : par_polar_box( thread ) ;
}

// return a double precision number
double
par_rng_dbl( const uint32_t thread )
{
#if (defined KISS_RNG)
  return par_KISS_dbl( thread ) ;
#elif (defined MWC_4096_RNG)
  return par_MWC_4096_dbl( thread ) ;
#elif (defined XOR_1024_RNG)
  return par_XOR_1024_dbl( thread ) ;
#elif (defined MWC_1038_RNG)
  return par_MWC_1038_dbl( thread ) ;
#else
  return par_WELL_512_dbl( thread ) ;
#endif
}

// accessor for ints
uint32_t
par_rng_int( const uint32_t thread )
{
  return (uint32_t)( UINT32_MAX * par_rng_dbl( thread ) ) ;
}

// read in the parallel rng state
int
read_par_rng_state( const char *infile )
{
  FILE *in = fopen( infile , "rb" ) ;
  struct QCDheader *get_header( FILE *__restrict in ) , * hdr ; 
  char *str ;

  if( in == NULL ) {
    fprintf( stderr , "[PAR_RNG] State file %s not found\n" , infile ) ;
    return GLU_FAILURE ;
  }
  
  if( ( hdr = get_header( in ) ) == NULL ) {
    fprintf( stderr , "[PAR_RNG] Unable to read state header\n" ) ;
    return GLU_FAILURE ;
  }

  // check Nthreads
  size_t Nthreads ;
  if( get_size_t( "NTHREADS" , hdr , &Nthreads ) == GLU_FAILURE ) {
    fprintf( stderr , "[PAR_RNG] NTHREADS not found in header\n" ) ;
    return GLU_FAILURE ;
  }
  if( Nthreads != (size_t)Latt.Nthreads ) {
    fprintf( stderr , "[PAR_RNG] RNG Nthreads not the same as Latt.Nthreads\n" ) ;
    return GLU_FAILURE ;
  }

  // figure out what RNG we are using and make sure it is consistent
  if( get_string( "RNG" , hdr , &str ) == GLU_FAILURE ) {
    fprintf( stderr , "[PAR_RNG] RNG type not found" ) ;
    return GLU_FAILURE ;
  }

#if (defined KISS_RNG)
  if( strcmp(  " PAR_KISS" , str ) ) {
    fprintf( stderr , "[PAR_RNG] state RNG differs from compiled (KISS) RNG\n" ) ;
    return GLU_FAILURE ;
  }
  if( read_par_KISS_table( in ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
#elif (defined MWC_4096_RNG)
  if( strcmp(  " PAR_MWC_4096" , str ) ) {
    fprintf( stderr , "[PAR_RNG] state RNG differs from compiled (MWC_4096) RNG\n" ) ;
    return GLU_FAILURE ;
  }
  if( read_par_MWC_4096_table( in ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
#elif (defined XOR_1024_RNG)
  if( strcmp(  " PAR_XOR_1024" , str ) ) {
    fprintf( stderr , "[PAR_RNG] state RNG differs from compiled (XOR_1024) RNG\n" ) ;
    return GLU_FAILURE ;
  }
  if( read_par_XOR_1024_table( in ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
#elif (defined MWC_1038_RNG)
  if( strcmp(  " PAR_MWC_1038" , str ) ) {
    fprintf( stderr , "[PAR_RNG] state RNG differs from compiled (MWC_1038) RNG\n" ) ;
    return GLU_FAILURE ;
  }
  if( read_par_MWC_1038_table( in ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
#else
  if( strcmp(  " PAR_WELL_512" , str ) ) {
    fprintf( stderr , "[PAR_RNG] state RNG differs from compiled (WELL_512) RNG\n" ) ;
    return GLU_FAILURE ;
  }
  if( read_par_WELL_512_table( in ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
#endif

  fclose( in ) ;

  return GLU_SUCCESS ;
}

// write out the parallel rng
void
write_par_rng_state( const char *outfile )
{
  FILE *out = fopen( outfile , "wb" ) ;

  // write out in a NERSC-like header so I can reuse that code
  fprintf( out , "BEGIN_HEADER\n" ) ;
  // write out number of threads
  fprintf( out , "NTHREADS = %u\n" , Latt.Nthreads ) ;
  // write out the rng type
#if (defined KISS_RNG)
  fprintf( out , "RNG = PAR_KISS\n" ) ;
#elif (defined MWC_4096_RNG)
  fprintf( out , "RNG = PAR_MWC_4096\n" ) ;
#elif (defined XOR_1024_RNG)
  fprintf( out , "RNG = PAR_XOR_1024\n" ) ;
#elif (defined MWC_1038_RNG)
  fprintf( out , "RNG = PAR_MWC_1038\n" ) ;
#else
  fprintf( out , "RNG = PAR_WELL_512\n" ) ;
#endif
  // tell us how big the table is
  fprintf( out, "RNG_TABLE = %d\n" , RNG_TABLE ) ;
  fprintf( out , "END_HEADER\n" ) ;

#if (defined KISS_RNG)
  write_par_KISS_table( out ) ;
#elif (defined MWC_4096_RNG)
  write_par_MWC_4096_table( out ) ;
#elif (defined XOR_1024_RNG)
  write_par_XOR_1024_table( out ) ;
#elif (defined MWC_1038_RNG)
  write_par_MWC_1038_table( out ) ;
#else
  write_par_WELL_512_table( out ) ;
#endif

  fclose( out ) ;

  return ;
}
