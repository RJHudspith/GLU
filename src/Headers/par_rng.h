/**
   @file par_rng.h
   @brief prototype functions for the parallel rng
 */
#ifndef PAR_RNG_H
#define PAR_RNG_H

void
free_par_rng( void ) ;

int
initialise_par_rng( const char *rng_file ) ;

void 
par_generate_NCxNC( GLU_complex U[ NCNC ] , 
		    const uint32_t thread ) ;

GLU_real
par_polar( const uint32_t thread ) ;

GLU_complex
par_polar_box( const uint32_t thread ) ;

double
par_rng_dbl( const uint32_t thread ) ;

uint32_t
par_rng_int( const uint32_t thread ) ;

int
read_par_rng_state( const char *infile ) ;

void
write_par_rng_state( const char *outfile ) ;

#endif
