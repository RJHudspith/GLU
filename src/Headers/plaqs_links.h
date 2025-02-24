/*
Copyright 2013-2025 Renwick James Hudspith

    This file (plaqs_links.h) is part of GLU.

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
   @file plaqs_links.h
   @brief function definitions for the gauge-fixing tests
 */
#ifndef GLU_PLAQS_LINKS_H
#define GLU_PLAQS_LINKS_H

/**
   @fn double all_plaquettes( const struct site *__restrict lat , double *__restrict sp_plaq , double *__restrict t_plaq ) 
   @brief computes the average plaquette for all directions
   @param lat :: gauge field
   @param sp_plaq :: average spatial plaquette
   @param t_plaq :: average temporal plaquette
   does not require fields to be SU(N) uses an identity for the specific trace of four matrices of the pattern \f$ U.U.U^{\dagger}.U^{\dagger}\f$

   @return 
   \f[

   \frac{1}{NC * ( ND - 1 ) * ( ND - 2 ) * LVOLUME} \sum_{x,\mu\neq\nu} \Re\left(tr\left[U_{\mu\nu}(x)\right]\right)

   \f]
 */
double 
all_plaquettes( const struct site *__restrict lat ,
		double *sp_plaq ,
		double *t_plaq ) ;

/**
   @fn double av_plaquette( const struct site *__restrict lat ) 
   @brief computes the average plaquette for all directions
   @param lat :: gauge field
   does not require fields to be SU(N) uses an identity for the specific trace of four matrices of the pattern \f$ U.U.U^{\dagger}.U^{\dagger}\f$

   @return 
   \f[

   \frac{1}{NC * ( ND - 1 ) * ( ND - 2 ) * LVOLUME} \sum_{x,\mu\neq\nu} \Re\left(tr\left[U_{\mu\nu}(x)\right]\right)

   \f]
 */
double 
av_plaquette( const struct site *__restrict lat ) ;

/**
   @fn double av_plaquette( const struct site *__restrict lat ) 
   @brief computes the average plaquette for all directions
   @param red :: reduction array. Sums into element 3 + #CLINE*Latt.Nthreads
   @param lat :: gauge field
   does not require fields to be SU(N) uses an identity for the specific trace of four matrices of the pattern \f$ U.U.U^{\dagger}.U^{\dagger}\f$
   @warning expected to be called in a parallel environment
 */
void 
av_plaquette_th( double *red ,
		 const struct site *__restrict lat ) ;

/**
   @fn double s_plaq( const struct site *__restrict lat ) 
   @brief computes the average plaquette in spatial direction 
   @param lat :: gauge field
   does not require fields to be SU(N)

   @return 
   \f[

   \frac{1}{NC * ( ND - 1 ) * ( ND - 2 ) * LVOLUME} \sum_{x,\mu\neq\nu,\mu,\nu \in ND-1} \Re\left(tr\left[U_{\mu\nu}(x)\right]\right)
   \f]
 */
double 
s_plaq( const struct site *__restrict lat ) ;

/**
   @fn double t_plaq( const struct site *__restrict lat ) 
   @brief computes the average plaquette in temporal direction 
   @param lat :: gauge field
   does not require fields to be SU(N)

   @return 
   \f[

   \frac{1}{NC * ( ND - 1 ) * LVOLUME} \sum_{x,\mu\neq (ND-1)} \Re\left(tr\left[ U_{\mu (ND-1)}(x)\right]\right)

   \f]
 */
double 
t_plaq( const struct site *__restrict lat ) ;

/**
   @fn int gauge_topological_meas( const struct site *__restrict lat , double *qtop_new , double *qtop_old , const int iter ) ;
   @brief computes the lattice topological charge
   @param red :: reduction array (GT in element 0 and Qtop in element 1)
   @param lat :: lattice fields
   @param qtop_new :: the new measurement of the topological charge
   @param qtop_old :: the old measurement of the topological charge
   @param iter :: iteration number
   @return #GLU_SUCCESS or #GLU_FAILURE as to whether we have met convergence criteria for finding an integer topological charge
*/
int
gauge_topological_meas_th( double *red ,
			   const struct site *__restrict lat ,
			   double *qtop_new ,
			   double *qtop_old ,
			   const size_t iter ) ;

/**
   @fn double all_links( const struct site *__restrict lat , double *__restrict sp_link , double *__restrict t_link ) 
   @brief computes the lattice average trace of all directions of links
   @param lat :: gauge field
   @param sp_link :: average spatial link trace
   @param t_link :: average temporal link trace
   @return
   \f[
   \frac{1}{ND*NC*LVOLUME} \sum_{x,\mu} \Re\left(tr\left[U_{\mu}(x)\right]\right)
   \f]
 */
double 
all_links( const struct site *__restrict lat ,
	   double *__restrict sp_link ,
	   double *__restrict t_link ) ;

/**
   @fn double indivlinks( const struct site *__restrict lat , GLU_real *max ) 
   @brief computes the lattice average trace of all directions of links
   @param lat :: gauge field
   @param max :: maximum link trace

   @return
   \f[
   \frac{1}{ND*NC*LVOLUME} \sum_{x,\mu} \Re\left(tr\left[U_{\mu}(x)\right]\right)
   \f]
 */
double
indivlinks( const struct site *__restrict lat , 
	    GLU_real *max ) ;

/**
   @fn double links( const struct site *__restrict lat ) 
   @brief computes the lattice average trace of all directions of links
   @param lat :: gauge field

   @return
   \f[
   \frac{1}{ND*NC*LVOLUME} \sum_{x,\mu} \Re\left(tr\left[U_{\mu}(x)\right]\right)
   \f]
 */
double 
links( const struct site *__restrict lat ) ;

/**
   @fn double s_links( const struct site *__restrict lat ) 
   @brief computes the lattice average trace of the spatial links
   @param lat :: gauge field

   @return
   \f[
   \frac{1}{(ND-1)*NC*LVOLUME} \sum_{x,mu\in ND-1} \Re\left(tr\left[U_{ND-1}(x)\right]\right)
   \f]
 */
double 
s_links( const struct site *__restrict lat ) ;

/**
   @fn double t_links( const struct site *__restrict lat ) 
   @brief computes the lattice average trace of the temporal links 
   @param lat :: gauge field
   @return
   \f[
   \frac{1}{NC*LVOLUME} \sum_{x} \Re\left(tr\left[U_{\mu}(x)\right]\right)
   \f]
 */
double
t_links( const struct site *__restrict lat ) ;

#endif
