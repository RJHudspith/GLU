#ifndef GLU_LINE_SEARCH_H
#define GLU_LINE_SEARCH_H

void
egauge_Landau( GLU_complex **gauge , 
	       const GLU_complex **in ,
	       const double alpha ) ;

void
exponentiate_gauge_CG( GLU_complex **gauge , 
		       const GLU_complex **in ,
		       const double alpha ) ;

/**
   @fn double approx_minimum( const size_t nmeas , const double alphas[ nmeas ] , const double functional[ nmeas ] )
   @brief finds the approximate minimum of alphas using GLU-bic splines
   @param nmeas :: number of measurements made
   @param alphas :: alphas tested
   @param functional :: the gauge functional at each alpha

   @return the alpha that approximately minimises the functional, or 0 
 */
double
approx_minimum( const size_t nmeas , 
		const double alphas[ nmeas ] ,
		const double functional[ nmeas ] ) ;

void
line_search_Coulomb( double *red ,
		     GLU_complex **gauge ,
		     const struct s_site *rotato ,
		     const struct draughtboard db ,
		     const struct site *lat ,
		     const GLU_complex **in ,
		     const size_t t ) ;

void
line_search_Landau( double *red ,
		    GLU_complex **gauge , 
		    const struct site *lat ,
		    const GLU_complex **in ) ;

#endif
