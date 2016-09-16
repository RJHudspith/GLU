/**
   @file relax.c
   @brief over relaxation code lives here
 */
#include "Mainfile.h"

#include "staples.h"    // all_staples()
#include "SU2_rotate.h" // rotation

// chroma's seems the cheapest at the moment
#define CHROMA_RELAX

// microcanonical su(2) update
void
microcanonical( GLU_complex *s0 ,
		GLU_complex *s1 )
{
  // chroma's is better as it has one fewer rotation
#ifdef CHROMA_RELAX
  // compute z
  register const double z = ( creal(*s0)*creal(*s0) - cimag(*s0)*cimag(*s0) -
			      creal(*s1)*creal(*s1) - cimag(*s1)*cimag(*s1) ) ;
  *s1 = -2. * creal( *s0 ) * ( *s1 ) ;
  *s0 = z - I * ( 2. * creal(*s0) * cimag(*s0) ) ;
#else
  // milc relax
  *s0 =  creal( *s0 ) - I * cimag( *s0 ) ;
  *s1 = -creal( *s1 ) - I * cimag( *s1 ) ;
#endif
  return ;
}

// overrelaxation algorithm
static void
overrelax( GLU_complex U[ NCNC ] , 
	   const GLU_complex staple[ NCNC ] )
{
  GLU_complex s0 , s1 ;
  double scale ;
  size_t i ;
  for( i = 0 ; i < NSU2SUBGROUPS ; i++ ) {
    only_subgroup( &s0 , &s1 , &scale , U , staple , i ) ;
    microcanonical( &s0 , &s1 ) ;
    #ifdef CHROMA_RELAX
    su2_rotate( U , s0 , s1 , i ) ;
    #else
    su2_rotate( U , s0 , s1 , i ) ;  
    su2_rotate( U , s0 , s1 , i ) ;
    #endif
  }
  return ;
}

// perform a heat-bath over the whole lattice
int
OR_lattice( struct site *lat ,
	    const struct draughtboard db )
{
  size_t i , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // update staples surrounding red links and update links
    #pragma omp parallel for private(i)
    for( i = 0 ; i < db.Nred ; i++ ) {
      GLU_complex stap[ NCNC ] GLUalign ;
      zero_mat( stap ) ;
      all_staples( stap , lat , db.red[i]  , mu , ND , SM_APE ) ;
      overrelax( lat[ db.red[i] ].O[mu] , stap ) ;
    }
    // update staples surrounding black links
    #pragma omp parallel for private(i)
    for( i = 0 ; i < db.Nblack ; i++ ) {
      GLU_complex stap[ NCNC ] GLUalign ;
      zero_mat( stap ) ;
      all_staples( stap , lat , db.black[i] , mu , ND , SM_APE ) ;
      overrelax( lat[ db.black[i] ].O[mu] , stap ) ;
    }
  }
  return GLU_SUCCESS ;
}

// clean up the over-relaxation
#ifdef CHROMA_RELAX
  #undef CHROMA_RELAX
#endif
