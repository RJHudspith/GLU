/**
   @file par_XOR_1024.h
   @brief prototype functions for the XOR_1024 rng
 */
#ifndef GLU_PAR_XOR_1024_H
#define GLU_PAR_XOR_1024_H

void
GLU_set_par_XOR_1024_table( const uint32_t seed[ Latt.Nthreads ] ) ;

void
free_par_XOR_1024( void ) ;

double
par_XOR_1024_dbl( const uint32_t thread ) ;

int
read_par_XOR_1024_table( FILE *rng_file ) ;

void
write_par_XOR_1024_table( FILE *rng_file ) ;

#endif
