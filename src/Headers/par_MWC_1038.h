#ifndef GLU_PAR_MWC_1038
#define GLU_PAR_MWC_1038

void
free_par_MWC_1038( void ) ;

void
GLU_set_par_MWC_1038_table( const uint32_t seed[ Latt.Nthreads ] ) ;

double 
par_MWC_1038_dbl( const uint32_t thread ) ;

int
read_par_MWC_1038_table( FILE *rng_file ) ;

void
write_par_MWC_1038_table( FILE *rng_file ) ;

#endif
