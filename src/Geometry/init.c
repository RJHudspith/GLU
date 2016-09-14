/**
   @file init.c
   @brief initialise the lattice
 */
#include "Mainfile.h"

#include "geometry.h"   // init_navig()
#include "SU2_rotate.h" // compute_pertinent_indices

// initialises the navigation for our lattice fields
void 
init_navig( struct site *__restrict lat )
{
  size_t i ; 
#pragma omp parallel for private(i)
  for(  i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu ;
    for(  mu = 0 ; mu < ND ; mu++ ) {
      lat[i].neighbor[mu] = gen_shift( i , mu ) ; 
      lat[i].back[mu] = gen_shift( i , -mu - 1 ) ;  
    }
  }
  return ;
}

// initialise the lattice information stored in Latt
void
init_latt( void )
{
  // these are neccessary for geometry and stuff ::  x,y,z,t geometry
  fprintf( stdout , "\n[DIMENSIONS] ( %zu x" , Latt.dims[0] ) ;
#if ND == 2
  Latt.Lcu = Latt.dims[0] ;
  Latt.Volume = Latt.Lcu * Latt.dims[1] ;
  fprintf( stdout , " %zu )\n\n" , Latt.dims[1] ) ; 
  return ;
#else
  size_t mu ;
  Latt.Lsq = Latt.dims[0] ; // defined :: LSQ
  for( mu = 1 ; mu < ND-2 ; mu++ ) {
    fprintf( stdout , " %zu x" , Latt.dims[ mu ] ) ;
    Latt.Lsq *= Latt.dims[mu] ;
  } 
  Latt.Lcu = Latt.Lsq * Latt.dims[ ND - 2 ] ;     // defined :: LCU
  Latt.Volume = Latt.Lcu * Latt.dims[ ND - 1 ] ;  // defined :: LVOLUME
  fprintf( stdout , " %zu x %zu )\n\n" , Latt.dims[mu] , Latt.dims[ mu+1 ] ) ; 
#endif
  // precompute the twiddle factors
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.twiddles[ mu ] = TWOPI / (double)Latt.dims[ mu ] ;
  }
  // set the global number of threads we will use
#if (defined _OPENMP ) && (defined HAVE_OMP_H)
  #pragma omp parallel
  { Latt.Nthreads = omp_get_num_threads( ) ; }
#else
  Latt.Nthreads = 1 ;
#endif
  fprintf( stdout , "[INIT] using %u thread(s) \n" , Latt.Nthreads ) ;
  // compute su2 indices
  fprintf( stdout , "[INIT] allocating su(2) subgroup indices\n" ) ;
  compute_pertinent_indices() ;
  return ;
}

// free anything allocated in the latt struct
void
free_latt( void )
{
  free( Latt.su2_data ) ;
}
