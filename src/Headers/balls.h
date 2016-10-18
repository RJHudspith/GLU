/**
   @file balls.c
   @brief GLU-balls calculation 
 */
#ifndef GLU_BALLS_H
#define GLU_BALLS_H

/**
   @fn int GLU_balls( struct site *lat , const struct sm_info SMINFO )
   @brief compute some glueballs
   @param lat :: lattice links
   @param SMINFO :: smearing information
 */
int
GLU_balls( struct site *lat ,
	   struct sm_info SMINFO ) ;

#endif
