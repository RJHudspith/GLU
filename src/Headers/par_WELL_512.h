/**
   @file par_WELL_512.h
   @brief prototype declarations for the WELL 512
 */
#ifndef GLU_PAR_WELL_512_H
#define GLU_PAR_WELL_512_H

void
free_par_WELL_512( void ) ;

void
GLU_set_par_WELL_512_table( const uint32_t seed[ Latt.Nthreads ] ) ;

double 
par_WELL_512_dbl( const uint32_t thread ) ;

int
read_par_WELL_512_table( FILE *rng_file ) ;

void
write_par_WELL_512_table( FILE *rng_file ) ;

#endif
