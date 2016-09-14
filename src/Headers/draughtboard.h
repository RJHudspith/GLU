/**
   @file draughtboard.h
   @brief prototype functions for setting the draughtboard
 */
#ifndef DRAUGHTBOARD_H
#define DRAUGHTBOARD_H

/**
   @fn void free_cb( struct draughtboard *db ) 
   @brief free the draughtboard
   @param db :: draughtboard structure
 */
void
free_cb( struct draughtboard *db ) ;

/**
   @fn void init_cb( struct draughtboard *db , const size_t LENGTH , const size_t DIR ) 
   @brief initialise the draughtboard
   @param db :: draughtboard structure
   @param LENGTH :: total length of the thing we are draughtboarding
   @param DIR :: number of directions in the geometry we are using
 */
void
init_cb( struct draughtboard *db ,
	 const size_t LENGTH ,
	 const size_t DIR ) ;

#endif
