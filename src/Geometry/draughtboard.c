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
