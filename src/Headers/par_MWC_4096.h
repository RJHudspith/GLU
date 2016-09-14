/**
   @file par_MWC_4096.h
   @brief prototype declarations for the MWC 4096
 */
#ifndef GLU_PAR_MWC_4096_H
#define GLU_PAR_MWC_4096_H

void
free_par_MWC_4096( void ) ;

void
GLU_set_par_MWC_4096_table( const uint32_t seed[ Latt.Nthreads ] ) ;

double 
par_MWC_4096_dbl( const uint32_t thread ) ;

int
read_par_MWC_4096_table( FILE *rng_file ) ;

void
write_par_MWC_4096_table( FILE *rng_file ) ;

#endif
