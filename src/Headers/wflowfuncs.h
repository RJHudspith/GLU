/*
    Copyright 2013 Renwick James Hudspith

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
   @var MEAS_START
   @brief at what point in the flow do we start measuring the topological stuff
**/
extern const double MEAS_START ;

/**
   @var TMEAS_STOP
   @brief at what point in the flow do we intend to stop
**/
extern const double TMEAS_STOP ;

/**
   @var WFLOW_STOP
   @brief what derivative value do we decide to (roughly) stop at
 */
extern const double WFLOW_STOP ;

/**
   @fn double deriv_euler( struct site *__restrict lat , double *flow , double *flow_next , const double t , const double delta_t )
   @brief Euler derivative
   @param lat :: lattice gauge field
   @param flow :: the flow before the new step
   @param flow_next :: the result from an integration step
   @param t :: flow time
   @param delta_t :: flow timestep
   @warning rewrites flow_next and flow
 */
double
deriv_euler( struct site *__restrict lat , 
	     double *flow , 
	     double *flow_next ,
	     const double t ,
	     const double delta_t ) ;

/**
   @fn double deriv_leapfrog( struct site *__restrict lat , double *flow_prev , double *flow , double *flow_next , const double t , const double delta_t )
   @brief computes the derivative using the leapfrog def
   @param lat :: lattice gauge field
   @param flow_prev :: the previous wilson flow result
   @param flow_next :: the result from an integration step
   @param flow :: the flow before the new step
   @param t :: flow time
   @param delta_t :: flow timestep
   @warning rewrites flow_prev, flow_next and flow
 */
double
deriv_leapfrog( struct site *__restrict lat , 
		double *flow_prev ,
		double *flow , 
		double *flow_next ,
		const double t ,
		const double delta_t ) ;

/**
   @fn void print_GG_info( const int SM_TYPE , const wflow_type WFLOW_TYPE )
   @brief prints out some relevant information for the flow routines
   @param SM_TYPE :: smearing type
   @param WFLOW_TYPE :: type of the wilson flow routine
 */
void
print_GG_info( const int SM_TYPE , 
	       const wflow_type WFLOW_TYPE ) ;

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
   @fn void step_distance( struct site *__restrict lat , struct spt_site *__restrict lat2 , struct spt_site_herm *__restrict Z , const double rk1 , const double rk2 , const double rk3 , const int SM_TYPE )
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
	       struct spt_site *__restrict lat2 ,
	       struct spt_site_herm *__restrict Z ,
	       const double rk1 ,
	       const double rk2 , 
	       const double rk3 , 
	       const int SM_TYPE ) ;

/**
   @fn void step_distance_memcheap( struct site *__restrict lat , struct spt_site *__restrict lat2 , struct spt_site *__restrict lat3 , struct spt_site *__restrict lat4 , struct spt_site_herm *__restrict Z , const double rk1 , const double rk2 , const double rk3 , const int SM_TYPE )
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
			struct spt_site *__restrict lat2 ,
			struct spt_site *__restrict lat3 ,
			struct spt_site *__restrict lat4 ,
			struct spt_site_herm *__restrict Z ,
			const double rk1 ,
			const double rk2 , 
			const double rk3 , 
			const int SM_TYPE ) ;

#endif
