/**
   @file relax.h
   @brief performs over-relaxation over the whole lattice
 */
#ifndef RELAX_H
#define RELAX_H

/**
   @fn void microcanonical( GLU_complex *s0 , GLU_complex *s1 )
   @brief microcanonically compute s0 and s1
   @param s0 :: (0,0) element of su2 matrix
   @param s1 :: (0,1) element of su2 matrix
 */
void
microcanonical( GLU_complex *s0 ,
		GLU_complex *s1 ) ;

/**
   @fn int OR_lattice( struct site *lat , const struct draughtboard db )
   @brief overrelax our links
 */
int
OR_lattice( struct site *lat ,
	    const struct draughtboard db ) ;

#endif
