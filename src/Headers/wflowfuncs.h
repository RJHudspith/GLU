/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (wflowfuncs.h) is part of GLU.

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
   @file wflowfuncs.h
   @brief prototype functions for the wilson flow routines
 */

#ifndef GLU_WFLOWFUNCS_H
#define GLU_WFLOWFUNCS_H

/**
   @var W0_STOP
   @brief what derivative value do we decide to (roughly) measure at
 */
extern const double W0_STOP ;

/**
   @var T0_STOP
   @brief what flow value do we decide to (roughly) perform measurements at
 */
extern const double T0_STOP ;

/**
   @fn const double evaluate_scale( double *der , const double *x , const double *meas , const size_t Nmeas , const double scale ,	const char *message )
   @brief evaluate the flow at scale via its cubic spline
   @param der :: high-order finite difference in GLU_splines.h
   @param x :: flow time
   @param meas :: flow measurement
   @param Nmeas :: number of flow measurements
   @param scale :: the scale the flow is evaluated at
   @param message :: identifier to make grepping easier
   @warning computes the derivative @der in this function
 */
double
evaluate_scale( double *der , 
		const double *x ,
		const double *meas ,
		const size_t Nmeas ,
		const double scale ,
		const char *message ) ;

/**
   @fn void print_flow( const struct wfmeas *curr , const double err , const double delta_t)
   @brief prints out the flow observable information
 */
void
print_flow( const struct wfmeas *curr ,
	    const double err ,
	    const double delta_t) ;

/**
   @fn void print_GG_info( void )
   @brief prints out some relevant information for the flow routines
 */
void
print_GG_info( void ) ;

/**
   @fn void scaleset( struct wfmeas *curr , const double W0 , const double T0 , const size_t count ) 
   @brief computes the flow in curr for T0 and W0
   @param curr :: measurement list
   @param T_0 :: wilson flow scale G(t) = W_0
   @param W_0 :: wilson flow scale t( dG(t)/dt ) = W_0
   @param count :: number of measurements made
 */
void
scaleset( struct wfmeas *curr , 
	  const double T_0 ,
	  const double W_0 ,
	  const size_t count ) ;

/**
   @fn void second_deriv( const double flow_prev , const double flow , const double flow_next , const double t , const double delta_t )
   @brief compute the dimensionless second derivative of the wilson flow
   @param flow_prev :: the previous wilson flow result
   @param flow_next :: the result from an integration step
   @param flow :: the flow before the new step
   @param t :: flow time
   @param delta_t :: flow timestep
   \f[
   t^3 \frac{\left( w( t + \delta t ) + w( t - \delta t ) - 2 w(t) \right)}{(\delta t)^2}
   \f]
 */
void
second_deriv( const double flow_prev ,
	      const double flow ,
	      const double flow_next ,
	      const double t ,
	      const double delta_t ) ;

/**
   @fn void step_distance( struct site *__restrict lat , struct s_site *__restrict lat2 , struct s_site_herm *__restrict Z , const double rk1 , const double rk2 , const double rk3 , const smearing_types SM_TYPE , void (*project)( GLU_complex log[ NCNC ] , GLU_complex *__restrict staple , const GLU_complex link[ NCNC ] , const double smear_alpha ) )
   @brief perform one rk4 wilson flow integration step
   @warning only SM_LOG and SM_STOUT available
   @param lat :: lattice gauge field
   @param lat2 :: temporary, the updated field for each RK4 step
   @param Z :: the antihermitian, lattice wide "Z" from the luescher paper
   @param rk1 :: the first rk4 stepping parameter
   @param rk2 :: the second rk4 stepping parameter
   @param rk3 :: the third rk4 stepping parameter
   @param SM_TYPE :: the smearing projection used 
   the parameters rk1 , rk2 and rk3 are as follows <br>
   rk1 = -9.0/17.0 * delta_t <br>
   rk2 = delta_t <br>
   rk3 = -delta_t <br>
 */
void
step_distance( struct site *__restrict lat ,
	       struct s_site *__restrict lat2 ,
	       struct s_site *__restrict Z ,
	       const double rk1 ,
	       const double rk2 , 
	       const double rk3 ,
	       const smearing_types SM_TYPE ,  
	       void (*project)( GLU_complex log[ NCNC ] , 
					 GLU_complex *__restrict staple , 
					 const GLU_complex link[ NCNC ] , 
					 const double smear_alpha ) ) ;

/**
   @fn void step_distance_memcheap( struct site *__restrict lat , struct s_site *__restrict lat2 , struct sp_site *__restrict lat3 , struct s_site *__restrict lat4 , struct s_site *__restrict Z , const double rk1 , const double rk2 , const double rk3 , const smearing_types SM_TYPE , void (*project)( GLU_complex log[ NCNC ] , GLU_complex *__restrict staple , const GLU_complex link[ NCNC ] , const double smear_alpha ) )
   @brief perform one rk4 wilson flow integration step
   @warning only SM_LOG and SM_STOUT available
   @param lat :: lattice gauge field
   @param lat2 :: temporary, the updated field on the previous slice
   @param lat3 :: temporary, the updated field on this slice
   @param lat4 :: temporary, the updated field on the very last slice
   @param Z :: the antihermitian, lattice wide "Z" from the luescher paper
   @param rk1 :: the first rk4 stepping parameter
   @param rk2 :: the second rk4 stepping parameter
   @param rk3 :: the third rk4 stepping parameter
   @param SM_TYPE :: the smearing projection used 
   the parameters rk1 , rk2 and rk3 are as follows <br>
   rk1 = -9.0/17.0 * delta_t <br>
   rk2 = delta_t <br>
   rk3 = -delta_t <br>
 */
void
step_distance_memcheap( struct site *__restrict lat ,
			struct s_site *__restrict lat2 ,
			struct s_site *__restrict lat3 ,
			struct s_site *__restrict lat4 ,
			struct s_site *__restrict Z ,
			const double rk1 ,
			const double rk2 , 
			const double rk3 , 
			const smearing_types SM_TYPE ,
			void (*project)( GLU_complex log[ NCNC ] , 
					 GLU_complex *__restrict staple , 
					 const GLU_complex link[ NCNC ] , 
					 const double smear_alpha ) ) ;

#endif
