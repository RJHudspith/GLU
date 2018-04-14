/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (clover.h) is part of GLU.

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
   @file clover.h 
   @brief Calculation of the naive field strength "G_{\mu\nu}" and topological charge
 */

#ifndef GLU_CLOVER_H
#define GLU_CLOVER_H

/**
   @fn void compute_Gmunu( double *__restrict GG , double *__restrict qtop , const struct site *__restrict lat ) ;
   @brief \f$ O(a^4) \f$ tree-improved field strength tensor from <a href="http://arxiv.org/abs/hep-lat/0203008"> paper </a>

   @param GG :: \f$ G_{\mu\nu}G_{\mu\nu} \f$
   @param qtop :: Naive topological charge
   @param lat :: lattice field
   @param mu :: direction mu
   @param nu :: orthogonal direction nu
   @param rho :: orthogonal direction again rho
   @param delta :: and the final (for 4D) orthogonal direction delta

   If we have not defined  <br>
   ANTIHERMITIAN ( traced antihermitian projection (BMW) ) <br>
   ANTIHERMITIAN_TRF ( traceless antihermitian projection (Chroma) ) <br>
   <br>
   we just use the traceless hermitian projection which is the same when <br>
   taking the products \f$ F_{\mu\nu}F_{\mu\nu} \f$ and  <br>
   \f$ e_{\mu\nu\rho\sigma}F_{\mu\nu}F_{\rho\sigma} \f$ <br>
   <br>
   As Jamie believes that the only way to incorporate the full \f$ O(a^4) \f$ <br>
   correction in the projector is to use the Log definition of the links <br>
   <br>
   #CLOVER_LOG_DEF ( traceless hermitian logarithm of the <br>
   exponentiated field strength tensor ) is also available.
 */
void
compute_Gmunu( double *__restrict GG , 
	       double *__restrict qtop ,
	       const struct site *__restrict lat ) ;

void
compute_Gmunu_th( double *red , 
		  const struct site *lat ) ;

/**
   @fn void compute_Gmunu_array( GLU_complex *__restrict *__restrict qtop , const struct site *__restrict lat )
   @brief \f$ O(a^4) \f$ tree-improved field strength tensor from <a href="http://arxiv.org/abs/hep-lat/0203008"> paper </a>
   @param qtop :: Naive topological charge matrix
   @param lat :: lattice field
 */
void
compute_Gmunu_array( GLU_complex *__restrict qtop ,
		     const struct site *__restrict lat ) ;

#endif
