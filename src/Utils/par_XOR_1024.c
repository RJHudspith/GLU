/**
   @file par_XOR_1024.c
   @brief parallel xor rng
 */
#include "Mainfile.h"

#include "GLU_bswap.h" // bswap64

static uint64_t **table ; // random number table
static uint32_t *xor_i ;  // xor 1024 index

// for the warm up
#define LCONG(a,i) ( 1812433253UL * ( a ^ ( a >> 30)) + i ) 

// returns a long long unsigned int
static uint64_t 
par_XOR_1024( const uint32_t thread ) 
{
  uint64_t s0 = table[ thread ][ xor_i[ thread ] ] ;
  uint64_t s1 = table[ thread ][ xor_i[ thread ] =			\
				 ( xor_i[ thread ] + 1 )&15 ] ;
  s1 ^= s1 << 31 ; // a
  s1 ^= s1 >> 11 ; // b
  s0 ^= s0 >> 30 ; // c
  return ( table[ thread ][ xor_i[ thread ] ] = s0 ^ s1 ) * 1181783497276652981LL ;
}

void
GLU_set_par_XOR_1024_table( const uint32_t seed[ Latt.Nthreads ] )
{
  table = malloc( Latt.Nthreads * sizeof( uint64_t* ) ) ;
  xor_i = malloc( Latt.Nthreads * sizeof( uint32_t ) ) ;
  size_t i , j ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {

    table[i] = malloc( RNG_TABLE * sizeof( uint64_t ) ) ;

    // type-pun the 64-bit table
    union {
      uint32_t theInts[2] ;
      uint64_t ulint ;
    } table64 ;

    table64.theInts[0] = seed[i] ;
    table64.theInts[1] = LCONG(seed[i],1) ;
    table[i][0] = table64.ulint ;

    // loop j
    for( j = 1 ; j < RNG_TABLE ; j++ ) {
      table64.theInts[0] = LCONG( table64.theInts[1] , 2*j     ) ;
      table64.theInts[1] = LCONG( table64.theInts[0] , 2*j + 1 ) ;
      table[i][j] = table64.ulint ;
    }

    xor_i[ i ] = RNG_TABLE - 1 ; // set p
  }
  return ;
}

void
free_par_XOR_1024( void )
{
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    free( table[i] ) ;
  }
  free( table ) ;
  free( xor_i ) ;
}

double
par_XOR_1024_dbl( const uint32_t thread ) 
{
  return ( par_XOR_1024( thread ) * 5.421010862427522e-20 ) ;
}

// read the MWC table
int
read_par_XOR_1024_table( FILE *rng_file )
{
  table = malloc( Latt.Nthreads * sizeof( uint64_t* ) ) ;
  xor_i = malloc( Latt.Nthreads * sizeof( uint32_t ) ) ;
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    table[i] = malloc( RNG_TABLE * sizeof( uint64_t ) ) ;
    if( fread( table[i] , sizeof( uint64_t ) , RNG_TABLE , rng_file ) != RNG_TABLE ) {
      fprintf( stderr , "[PAR_RNG] file read failure\n" ) ;
      return GLU_FAILURE ;
    }
  }
  if( fread( xor_i , sizeof( uint32_t ) , Latt.Nthreads , rng_file ) != Latt.Nthreads ) {
    fprintf( stderr , "[PAR_RNG] file read failure\n" ) ;
    return GLU_FAILURE ;
  }
  // swap them if we aren't big endian
  if( !WORDS_BIGENDIAN ) { 
    for( i = 0 ; i < Latt.Nthreads ; i++ ) {
      bswap_64( RNG_TABLE , table[i] ) ;
    }
    bswap_32( Latt.Nthreads , xor_i ) ;
  }
  return GLU_SUCCESS ;
}

// write the table and periferies
void
write_par_XOR_1024_table( FILE *rng_file )
{
  // write out big endian binary data
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    // put into big endian
    if( !WORDS_BIGENDIAN ) { 
      bswap_64( RNG_TABLE , table[i] ) ; 
    }
    // write the table
    fwrite( table[i] , sizeof( uint64_t ) , RNG_TABLE , rng_file ) ;
    // reswap
    if( !WORDS_BIGENDIAN ) { 
      bswap_64( RNG_TABLE , table[i] ) ; 
    }
  }
  // write out the carries
  if( !WORDS_BIGENDIAN ) { 
    bswap_32( Latt.Nthreads , xor_i ) ;
  }
  fwrite( xor_i , sizeof( uint32_t ) , Latt.Nthreads , rng_file ) ;
  if( !WORDS_BIGENDIAN ) { 
    bswap_32( Latt.Nthreads , xor_i ) ;
  }
  return ;
}
