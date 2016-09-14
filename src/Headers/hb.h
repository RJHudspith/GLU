/**
   @file hb.h
   @brief prototype declarations for heatbath algorithm
 */
#ifndef HB_H
#define HB_H

/**
   @fn int hb_lattice( struct site *lat , struct site *staple , const double invbeta , const struct draughtboard db )
   @brief perform heat-bath updates
 */
int
hb_lattice( struct site *lat ,
	    struct site *staple ,
	    const double invbeta ,
	    const struct draughtboard db ) ;

#endif
