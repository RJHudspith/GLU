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
  free( db -> red ) ;
  free( db -> black ) ;
}

// initialise the draughtboarding
// warning :: red and black are allocated in here!
void
init_cb( struct draughtboard *db ,
	 const size_t LENGTH ,
	 const size_t DIR ) 
{
  size_t i , nred = 0 , nblack = 0 ;
  int n[ ND ] ;
  for( i = 0 ; i < LENGTH ; i++ ) {
    const int mode_sum = get_mom_2piBZ( n , i , DIR ) ;
    if( mode_sum&1 ) {
      nred++ ;
    } else {
      nblack++ ;
    }
  }
  // malloc and set
  db -> red   = malloc( nred  * sizeof( size_t ) ) ;
  db -> black = malloc( nblack * sizeof( size_t ) ) ;
  nred = nblack = 0 ;
  // set back to zero
  for( i = 0 ; i < LENGTH ; i++ ) {
    const int mode_sum = get_mom_2piBZ( n , i , DIR ) ;
    if( mode_sum&1 ) {
      db -> red[ nred ] = i ;
      nred++ ;
    } else {
      db -> black[ nblack ] = i ;
      nblack++ ;
    }
  }
  db -> Nred = nred ; 
  db -> Nblack = nblack ;
  printf( "\n[DRAUGHTBOARD] initialised\n\n" ) ;

  return ;
}
