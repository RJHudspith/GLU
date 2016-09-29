/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (draughtboard.c) is part of GLU.

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
   @file draught_board.c
   @brief draught-boarding routines
 */
#include "Mainfile.h"
#include "geometry.h"

// free the draughtboard
void
free_cb( struct draughtboard *db )
{
  size_t i ;
  for( i = 0 ; i < db -> Ncolors ; i++ ) {
    free( db -> square[i] ) ;
  }
  free( db -> Nsquare ) ;
  free( db -> square ) ;
  return ;
}

// initialise the draughtboarding
// warning :: red and black are allocated in here!
// idea :: checkerboard the even sites and update the odd
int
init_cb( struct draughtboard *db ,
	 const size_t LENGTH ,
	 const size_t DIR ) 
{
  // counters and such
  size_t i ;
  int n[ ND ] ;
  size_t sum = 0 ;
  for( i = 0 ; i < ND ; i++ ) {
    if( ( Latt.dims[i]&1 ) == 1 ) {
      sum++ ;
    }
  }
  // complain
  if( sum != 0 && sum != ND ) {
    fprintf( stderr , "[DRAUGHTBOARD] can only generate all odd or all even\n" ) ;
    return GLU_FAILURE ;
  }

#ifdef IMPROVED_SMEARING
  db -> Ncolors = 32 ;
#else
  db -> Ncolors = ( sum == ND ) ? 3 : 2 ;
#endif
  db -> square = malloc( db -> Ncolors * sizeof( size_t* ) ) ;
  db -> Nsquare = malloc( db -> Ncolors * sizeof( size_t ) ) ;

  // set the draughtboard numbers to zero
  for( i = 0 ; i < db -> Ncolors ; i++ ) {
    db -> Nsquare[i] = 0 ;
  }

  // compute how many sites for each color, I understand they will be 
  // LVOLUME/db->Ncolors but we need to be careful about the exact numbers for
  // each bin when doing IMPROVED_SMEARING
  for( i = 0 ; i < LENGTH ; i++ ) {
    const size_t midx = ((size_t)get_mom_2piBZ( n , i , DIR ))%(db -> Ncolors) ;
    db -> Nsquare[ midx ]++ ;
  }

  // allocate the coloring and set the Nsquares to zero
  for( i = 0 ; i < db -> Ncolors ; i++ ) {
    //#ifdef verbose
    printf( "[DRAUGHTBOARD] COLOR_%zu :: %zu \n" , i , db -> Nsquare[i] ) ;
    //#endif
    db -> square[i] = malloc( db -> Nsquare[i] * sizeof( size_t ) ) ;
    db -> Nsquare[i] = 0 ;
  }

  // set back to zero and redo recording the index of each
  for( i = 0 ; i < LENGTH ; i++ ) {
    const size_t midx = ((size_t)get_mom_2piBZ( n , i , DIR ))%(db -> Ncolors) ;
    db -> square[ midx ][ db -> Nsquare[ midx ] ] = i ;
    db -> Nsquare[ midx ]++ ;
  }

  printf( "\n[DRAUGHTBOARD] initialised\n\n" ) ;

  return GLU_SUCCESS ;
}

int
test_db( struct site *lat ,
	 const struct draughtboard db )
{
  // test that each point on the db has a neighbour
  size_t i , j , c , mu ;
  for( i = 0 ; i < db.Nsquare[0] ; i++ ) {
    for( c = 1 ; c < db.Ncolors ; c++ ) {
      for( j = 0 ; j < db.Nsquare[c] ; j++ ) {
	if( db.square[0][i] == db.square[c][j] ) {
	  printf( "Fucked \n" ) ;
	  return GLU_FAILURE ;
	}
      }
    }
  }
  return GLU_SUCCESS ;
}
