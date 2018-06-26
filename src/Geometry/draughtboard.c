/*
    Copyright 2013-2018 Renwick James Hudspith

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

#include "geometry.h" //   get_mom_2piBZ()

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

// get the correct index
static size_t
get_midx( const size_t i ,
	  const size_t DIR ,
	  const size_t Ncolors )
{
  int n[ ND ] ;
  get_mom_2piBZ( n , i , DIR ) ;
  register size_t even_sum = 0 , odd_sum = 0 , mu ;
  for( mu = 0 ; mu < DIR ; mu++ ) {
    if( ( Latt.dims[mu]&1 ) == 0 ) {
      even_sum += n[mu] ;
    } else {
      odd_sum  += n[mu] ;
    }
  }
#ifdef verbose
  fprintf( stdout , "(%zu) even :: %zu || odd %zu -> %zu \n" , 
	   i , even_sum , odd_sum , ( even_sum%2 + odd_sum )%Ncolors ) ;
#endif
  return ( even_sum%2 + odd_sum )%Ncolors ;
}

// initialise the draughtboarding
// warning :: red and black are allocated in here!
// idea :: checkerboard the even sites and update the odd
int
init_cb( struct draughtboard *db ,
	 const size_t LENGTH ,
	 const size_t DIR ) 
{
  size_t i ;

  // initially set up the database
  db -> Ncolors = 3 ;
  db -> square  = malloc( db -> Ncolors * sizeof( size_t* ) ) ;
  db -> Nsquare = malloc( db -> Ncolors * sizeof( size_t ) ) ;
  
  // set the counters to zero
  for( i = 0 ; i < db -> Ncolors ; i++ ) {
    db -> Nsquare[i] = 0 ;
  }

  // get the coloring
  for( i = 0 ; i < LENGTH ; i++ ) {
    db -> Nsquare[ get_midx( i , DIR , db -> Ncolors ) ]++ ;
  }

  // if we don't have an odd dim we shorten this
  if( ( db -> Nsquare[ 2 ] ) == 0 ) {
    db -> Ncolors = 2 ;
  }

  // allocate the coloring and set the Nsquares to zero
  for( i = 0 ; i < db -> Ncolors ; i++ ) {
    #ifdef verbose
    fprintf( stdout , "[DB] COLOR_%zu :: %zu \n" , i , db -> Nsquare[i] ) ;
    #endif
    db -> square[i] = malloc( db -> Nsquare[i] * sizeof( size_t ) ) ;
    db -> Nsquare[i] = 0 ;
  }

  // set back to zero and redo recording the index of each
  for( i = 0 ; i < LENGTH ; i++ ) {
    const size_t midx = get_midx( i , DIR , db -> Ncolors ) ;
    db -> square[ midx ][ db -> Nsquare[ midx ] ] = i ;
    db -> Nsquare[ midx ]++ ;
  }

  fprintf( stdout , "\n[DB] initialised\n\n" ) ;

  return GLU_SUCCESS ;
}

int
test_db( struct site *lat ,
	 const struct draughtboard db )
{
  // test that each point on the db has a neighbour
  size_t i , j , c ;
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
