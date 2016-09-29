/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (init.c) is part of GLU.

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
  size_t mu ;
  fprintf( stdout , "\n[DIMENSIONS] ( %zu x" , Latt.dims[0] ) ;
#if ND == 2
  Latt.Lcu = Latt.dims[0] ;
  Latt.Volume = Latt.Lcu * Latt.dims[1] ;
  fprintf( stdout , " %zu )\n\n" , Latt.dims[1] ) ; 
#else
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
