/**
   @file SU2_rotate.h
   @brief prototype declarations for su2 hits
 */
#ifndef SU2_ROTATE_H
#define SU2_ROTATE_H

/**
   @fn void compute_pertinent_indices( void )
   @brief computes the various su2 subgoup indices of our NCxNC matrix
   @warning allocates memory

   I would like to thank Chroma, and in particular Edwards for this code.
 */
void
compute_pertinent_indices( void ) ;

/**
   @fn void free_su2_data( void )
   @brief frees the su2 subgroup struct allocated in givens.c
 */
void
free_su2_data( void ) ;

/**
   @fn void only_subgroup( GLU_complex *s0 , GLU_complex *s1 , double *scale , const GLU_complex U[ NCNC ] , const GLU_complex staple[ NCNC ] , const size_t su2_index )
   @brief gives the action for a specific subgroup efficiently by only computing the contribution for the specific index
   @param s0 :: first element of the su2 matrix
   @param s1 :: second element of the su2 matrix
   @param scale :: inverse determinant of the su2 matrix
   @param U :: link matrix
   @param staple :: surrounding staples of link @U
   @param su2_index :: su2 subgroup index
*/
void
only_subgroup( GLU_complex *s0 ,
	       GLU_complex *s1 ,
	       double *scale ,
	       const GLU_complex U[ NCNC ] ,
	       const GLU_complex staple[ NCNC ] ,
	       const size_t su2_index ) ;

/**
   @fn void shortened_su2_multiply( GLU_complex *w , const GLU_complex a , const GLU_complex b , const GLU_complex c , const GLU_complex d , const size_t su2_index )
   @brief su(2) multiply w = su2[ su2_index ] * w
   @param w :: matrix being hit on the left
   @param a :: top left su(2) element
   @param b :: top right su(2) element
   @param c :: bottom left su(2) element
   @param d :: bottom right su(2) element
   @param su2_index :: index describing which subgroup we are using
 */
void
shortened_su2_multiply( GLU_complex *w , 
			const GLU_complex a , 
			const GLU_complex b , 
			const GLU_complex c , 
			const GLU_complex d , 
			const size_t su2_index ) ;

/**
   @fn void shortened_su2_multiply_dag( GLU_complex *U , const GLU_complex a , const GLU_complex b , const GLU_complex c , const GLU_complex d , const size_t su2_index )
   @brief su(2) multiply U = U * su2[ su2_index ] ^{dagger}
   @param w :: matrix being hit on the right by a daggered SU(2)
   @param a :: top left su(2) element
   @param b :: top right su(2) element
   @param c :: bottom left su(2) element
   @param d :: bottom right su(2) element
   @param su2_index :: index describing which subgroup we are using
 */
void
shortened_su2_multiply_dag( GLU_complex *U , 
			    const GLU_complex a , 
			    const GLU_complex b , 
			    const GLU_complex c , 
			    const GLU_complex d , 
			    const size_t su2_index ) ;

/**
   @fn void su2_rotate( GLU_complex U[ NCNC ] , const GLU_complex s0 , const GLU_complex s1 , const size_t su2_index ) ;
   @brief rotate a link U by an su2 matrix
   @param U :: link matrix
   @param s0 :: first element of su(2) matrix
   @param s1 :: second element of su(2) matrix
   @param su2_index :: subgroup index being used in the rotation

   Computes \f$ U \rightarrow s_{\text{su2_index}} U\f$.
 */
void
su2_rotate( GLU_complex U[ NCNC ] ,
	    const GLU_complex s0 ,
	    const GLU_complex s1 ,
	    const size_t su2_index ) ;

#endif
