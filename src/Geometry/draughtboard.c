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

#include <assert.h>

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
    if( db -> Nsquare[i] == 0 ) {
      db -> square[i] = NULL ;
      return GLU_FAILURE ;
    } else {
      db -> square[i] = malloc( db -> Nsquare[i] * sizeof( size_t ) ) ;
    }
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

// very different version required for rectangle actions
int
init_improved_cb( struct draughtboard *db )
{
  assert( ND == 4 ) ;
  for( int mu = 0 ; mu < ND ; mu++ ) {
    assert( Latt.dims[mu]%4 == 0 && Latt.dims[mu] >=4 ) ;
  }
#ifdef SYMANZIK_ONE_LOOP
  assert( true ) ;
#endif
  
  size_t i ;

  // initially set up the database
  db -> Ncolors = 16 ;
  db -> square  = malloc( db -> Ncolors * sizeof( size_t* ) ) ;
  db -> Nsquare = malloc( db -> Ncolors * sizeof( size_t ) ) ;
  
  // set the counters to zero
  for( i = 0 ; i < db -> Ncolors ; i++ ) {
    db -> Nsquare[i] = LVOLUME/4 ;
    db->square[ i ] = malloc( db->Nsquare[i]*sizeof( size_t ) ) ;
    memset( db->square[i] , 0 , sizeof( size_t )*db->Nsquare[i] ) ;
  }

  const size_t multipliers[ ND ] = { Latt.dims[0]/4 , Latt.dims[1]/4 , Latt.dims[2]/4 , Latt.dims[3]/4 } ;
  
  // get the coloring
  for( size_t mu = 0 ; mu < ND ; mu++ ) {

    // first orthogonal direction is either x or y when mu == x -> arbitrary choice
    const size_t orth = ( mu == 0 ) ? 1 : 0 ;

    // ones we can use are a knights-move away
    int knightmove[ ND ] = {} ;
    knightmove[mu] = 1 ; knightmove[ orth ] = 2 ;

    // in the mu direction we only need to move 2 away
    int shift_mu[ ND ] = {} ;
    shift_mu[mu] = 2 ;

    // simple 4x4 pattern from origin mu == 0, and in arbitrary ortho direction that we choose to be y
    size_t leaf[ 4 ][ 4*multipliers[orth]*multipliers[mu] ] ;

    for( size_t phase = 0 ; phase < 4 ; phase++ ) {
      int off[4] = {} ;
      switch( phase ) {
      case 0 :
	off[ mu ] = 0 ; off[ orth ] = 0 ;
	break ;
      case 1 :
	off[ mu ] = 1 ; off[ orth ] = 1 ;
	break ;
      case 2 :
	off[ mu ] = 1 ; off[ orth ] = 0 ;
	break ;
      case 3 :
	off[ mu ] = 0 ; off[ orth ] = 1 ;
	break ;
      }
      fprintf( stdout , "mu == %zu || leaf[%zu] = ( " , mu , phase ) ;
      size_t counter = 0 ;
      for( size_t outer = 0 ; outer < 2*multipliers[orth] ; outer++ ) {
	if( outer == 0 ) {
	  leaf[ phase ][ counter ] = compute_spacing( off , 0 , ND )  ; // defines our starting point
	} else {
	  leaf[ phase ][ counter ] = compute_spacing( knightmove , leaf[ phase ][counter-1] , ND ) ;
	}
	fprintf( stdout , "%zu " , leaf[phase][counter] ) ;
	counter++ ;
	for( size_t inner = 1 ; inner < 2*multipliers[mu] ; inner++  ) {
	  leaf[ phase ][ counter ] = compute_spacing( shift_mu , leaf[phase][counter-1] , ND )  ; // defines our starting point
	  fprintf( stdout , "%zu " , leaf[phase][counter] ) ;
	  counter++ ;
	}
      }
      fprintf( stdout , ")\n" ) ;
      fprintf( stdout , "Counter ---> %zu\n" , counter ) ;
    }
    size_t d[2] = {} , counter = 0;
    for( size_t nu = 1 ; nu < ND ; nu++ ) {
      if( nu != mu && nu != orth ) {
	d[counter++] = nu ; 
      }
    }
    fprintf( stdout , "Ortho dirs == %zu %zu %zu\n\n" , orth , d[0] , d[1] ) ;
    
    for( size_t phase = 0 ; phase < 4 ; phase++ ) {
      const size_t idx = phase + 4*mu ;
      // loop two similar directions
      size_t i = 0 ;
      for( size_t orth1 = 0 ; orth1 < 4*multipliers[d[1]] ; orth1++ ) {
	for( size_t orth2 = 0 ; orth2 < 4*multipliers[d[0]] ; orth2++ ) {
	  int sep[4] = { } ;
	  sep[ d[0] ] = orth2 ; sep[ d[1] ] = orth1 ;
	  const size_t offset = compute_spacing( sep , 0 , ND ) ;
	  for( size_t munu = 0 ; munu < 4*multipliers[orth]*multipliers[mu] ; munu++ ) {
	    db->square[idx][i] = leaf[ (orth1+orth2+phase)%4 ][ munu ] + offset ;
	    i++ ;
	  }
	}
      }
    }
  }
  fprintf( stdout , "\n[DB] initialised\n\n" ) ;
  return GLU_SUCCESS ;
}
