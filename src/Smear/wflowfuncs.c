/*
Copyright 2013-2025 Renwick James Hudspith

    This file (wflowfuncs.c) is part of GLU.

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
   @file wflowfuncs.c
   @brief helper functions for the wilson flow
 */
#include "Mainfile.h"

#include "clover.h"      // compute_Gmumu_th
#include "GLU_splines.h" // spline evaluations
#include "plaqs_links.h" // clover discretisation
#include "projectors.h"  // stout projection
#include "staples.h"     // computes staples

// controls for the wilson flow these get externed!
const double W0_STOP = 0.1125*(NCNC-1)/(double)NC ; // BMW's choice for the W_0 parameter
const double T0_STOP = 0.1125*(NCNC-1)/(double)NC ; // Martin Luescher's choice for the t_0 scale

// shortening function needs to be out of alphabetical order because
// it is called by flow directions
static inline void
make_short_log( GLU_complex short_staple[ HERMSIZE ] , 
		const GLU_complex staple[ NCNC ] )
{
#if NC == 3
  short_staple[0] = staple[ 0 ] ;
  short_staple[1] = staple[ 1 ] ;
  short_staple[2] = staple[ 2 ] ;
  short_staple[3] = staple[ 4 ] ;
  short_staple[4] = staple[ 5 ] ;
#elif NC == 2
  short_staple[0] = staple[ 0 ] ;
  short_staple[1] = staple[ 1 ] ;
#else
  size_t i, j , idx = 0 ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = i ; j < NC ; j++ ) {
      short_staple[ idx ] = staple[ j + i*NC ] ;
      idx ++ ;
    }
  }
#endif
  return ;
}

// accumulate the Z-matrix 
static inline void
set_zmatrix( struct s_site *__restrict Z ,
	     const GLU_complex short_staple[ HERMSIZE ] ,
	     const double multiplier ,
	     const size_t i ,
	     const size_t mu ) 
{
#if NC == 3
  Z[i].O[mu][0] += multiplier * creal( short_staple[0] ) +\
    I * multiplier * creal( short_staple[3] ) ;
  Z[i].O[mu][1] += multiplier * short_staple[1] ;
  Z[i].O[mu][2] += multiplier * short_staple[2] ;
  Z[i].O[mu][3] += multiplier * short_staple[4] ;
#elif NC == 2
  Z[i].O[mu][0] += multiplier * creal( short_staple[0] ) ;
  Z[i].O[mu][1] += multiplier * short_staple[1] ;
#else
  size_t elements ;
  for( elements = 0 ; elements < TRUE_HERM ; elements++ ) {
    Z[i].O[mu][ elements ] += multiplier * short_staple[elements] ;
  }
#endif
}

// wilson flow in all ND directions ...
static void
flow_directions( struct s_site *__restrict lat2 ,
		 struct s_site *__restrict Z ,
		 const struct site *__restrict lat ,
		 const double multiplier ,
		 const double delta_t ,
		 const size_t i ,
		 const size_t it ,
		 const size_t mu ,
		 const smearing_types SM_TYPE ,
		 void (*project)( GLU_complex log[ NCNC ] , 
				  GLU_complex *__restrict staple , 
				  const GLU_complex link[ NCNC ] , 
				  const double smear_alpha ) )
{
  GLU_complex staple[ NCNC ] GLUalign = {} , temp[ NCNC ] GLUalign ;
  GLU_complex short_staple[ HERMSIZE ] GLUalign ; // does not need to be inited
  // first element
  #ifdef IMPROVED_SMEARING
  all_staples_improve( staple , lat , it , mu , ND , SM_TYPE ) ;
  #else
  all_staples( staple , lat , it , mu , ND , SM_TYPE ) ;
  #endif
  // default to STOUT ...
  switch( SM_TYPE ) {
  case SM_LOG :
    make_short_log( short_staple , staple ) ;
    break ;
  default :
    multab_dag( temp , staple , lat[it].O[mu] ) ; 
    Hermitian_proj_short( short_staple , temp ) ;
    break ;
  }
  // poke back into the accumulated Z-matrix
  set_zmatrix( Z , short_staple , multiplier , it , mu ) ;
  // perform the stout projection on the hermitian-shortened links
  project( lat2[i].O[mu] , Z[it].O[mu] , lat[it].O[mu] , delta_t ) ; 
  return ;
}

// memory-expensive runge-kutta step
static void
RK3step( struct wflow_temps WF ,
	 struct site *__restrict lat ,
	 const double multiplier , 
	 const double delta_t ,
	 const smearing_types SM_TYPE ,
	 void (*project)( GLU_complex log[ NCNC ] , 
			  GLU_complex *__restrict staple , 
			  const GLU_complex link[ NCNC ] , 
			  const double smear_alpha ) )
{
  size_t i , mu ;
#pragma omp for private(i) collapse(2) SCHED 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      GLU_complex staple[ NCNC ] GLUalign = {} , temp[ NCNC ] GLUalign = {} ;
      GLU_complex short_staple[ HERMSIZE ] GLUalign = {} ;
      switch( SM_TYPE ) {
      case SM_LOG :
        #ifdef IMPROVED_SMEARING
	all_staples_improve( staple , lat , i , mu , ND , SM_TYPE ) ;
        #else
	all_staples( staple , lat , i , mu , ND , SM_TYPE ) ;
        #endif
	// log smearing just requires an set the log
	make_short_log( short_staple , staple ) ;
	set_zmatrix( WF.Z , short_staple , multiplier , i , mu ) ;
	project( WF.lat2[i].O[mu] , WF.Z[i].O[mu] , lat[i].O[mu] , delta_t ) ; 
	break ;
      default :
        #ifdef IMPROVED_SMEARING
	all_staples_improve( staple , lat , i , mu , ND , SM_TYPE ) ;
        #else
	all_staples( staple , lat , i , mu , ND , SM_TYPE ) ;
        #endif
	// hermitian projection
	multab_dag( temp , staple , lat[i].O[mu] ) ; 
	Hermitian_proj_short( short_staple , temp ) ;
	set_zmatrix( WF.Z , short_staple , multiplier , i , mu ) ;
	project( WF.lat2[i].O[mu] , WF.Z[i].O[mu] , lat[i].O[mu] , delta_t ) ;
	break ;
      }
    }
  }
  // copy into lat
#pragma omp for private(i) collapse(2)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      equiv( lat[i].O[mu] , WF.lat2[i].O[mu] ) ;
    }
  }
  return ;
}

// computes one of the RK steps, doesn't matter which one
static void
RK3step_memcheap( struct wflow_temps WF ,
		  struct site *__restrict lat ,
		  const double multiplier ,
		  const double step ,
		  const smearing_types SM_TYPE ,
		  void (*project)( GLU_complex log[ NCNC ] , 
				   GLU_complex *__restrict staple , 
				   const GLU_complex link[ NCNC ] , 
				   const double smear_alpha ) )
{
  size_t i ;
#ifdef IMPROVED_SMEARING
  const size_t back = lat[ lat[0].back[ ND-1 ] ].back[ ND-1 ] ;
#else
  const size_t back = lat[0].back[ ND - 1 ] ;
#endif
  // split volume - wise
  #pragma omp for private(i) SCHED
  #ifdef IMPROVED_SMEARING
  for( i = 0 ; i < 2*LCU ; i++ ) {
  #else
  for( i = 0 ; i < LCU ; i++ ) {
  #endif
    const size_t bck = back + i ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      flow_directions( WF.lat4 , WF.Z , lat , multiplier , 
		       step , i , bck , mu , SM_TYPE , project ) ; 
    }
  }
  size_t t ;
  #ifdef IMPROVED_SMEARING
  for( t = 0 ; t < Latt.dims[ ND - 1 ] - 2 ; t++ ) {
  #else
  for( t = 0 ; t < Latt.dims[ ND - 1 ] - 1 ; t++ ) {
  #endif
    const size_t slice = LCU * t ;
    // perform the wilson flow for a point on the slice
    #pragma omp for private(i) SCHED
    for( i = 0 ; i < LCU ; i++ ) {
      const size_t it = slice + i ;
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	flow_directions( WF.lat2 , WF.Z , lat , multiplier , 
			 step , i , it , mu , SM_TYPE , project ) ;
      }
    }
    // swap over the temporary lattice fields
#ifdef IMPROVED_SMEARING
    const size_t bck = lat[ lat[ slice ].back[ ND-1 ] ].back[ ND-1 ] ;
#else
    const size_t bck = lat[ slice ].back[ ND -1 ] ;
#endif
#pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      //put temp into the previous time-slice
      register size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
        #ifdef IMPROVED_SMEARING
	if( t > 1 ) { 
	  register const size_t back = bck + i ;
	  equiv( lat[back].O[mu] , WF.lat3[i].O[mu] ) ;
	}
	// put the evaluation two time slices before into the front half of lat3
	// and the evaluation one time slice before into the lower half
	equiv( WF.lat3[i].O[mu] , WF.lat3[i+LCU].O[mu] ) ;
	equiv( WF.lat3[i+LCU].O[mu] , WF.lat2[i].O[mu] ) ;
        #else
	if( t != 0 ) { 
	  register const size_t back = bck + i ;
	  equiv( lat[back].O[mu] , WF.lat3[i].O[mu] ) ;
	}
	//make temporary lat3 lat2 again and repeat
	equiv( WF.lat3[i].O[mu] , WF.lat2[i].O[mu] ) ;
        #endif
      }
      //
    }
  }
  // put the last couple back in ....
  const size_t slice = LCU * t ;
  const size_t behind = lat[ slice ].back[ ND - 1 ] ;
#ifdef IMPROVED_SMEARING
  const size_t behind2 = lat[ behind ].back[ ND-1 ] ;
#endif
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    register const size_t back = behind + i ;
    register const size_t it = slice + i ; 
    register size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      #ifdef IMPROVED_SMEARING
      equiv( lat[behind2+i].O[mu] , WF.lat3[i].O[mu] ) ; 
      equiv( lat[back].O[mu]      , WF.lat3[i+LCU].O[mu] ) ; 
      equiv( lat[it].O[mu]        , WF.lat4[i].O[mu] ) ; 
      equiv( lat[it+LCU].O[mu]    , WF.lat4[i+LCU].O[mu] ) ; 
      #else
      equiv( lat[back].O[mu] , WF.lat3[i].O[mu] ) ; 
      equiv( lat[it].O[mu]   , WF.lat4[i].O[mu] ) ; 
      #endif
    }
  }
  return ;
}

// evaluate the flow using splines, evaluate a a particular scale
static double
evaluate_scale( double *der , 
		const double *x ,
		const double *meas ,
		const size_t Nmeas ,
		const double scale )
{
  // compute spline for observable meas, x must be sorted smallest to largest
  size_t i , change_up = 0 ;
  for( i = 0 ; i < Nmeas ; i++ ) {
    spline_derivative( der , x , meas , Nmeas ) ;
    // W_0 scale
    if( i > 0 && ( meas[ i ] > scale ) && ( meas[ i - 1 ] < scale ) ) {
      change_up = i ;
    }
    #ifdef verbose
    fprintf( stdout , "%g %g \n" , x[ i ] , meas[ i ] ) ;
    #endif
  }
  #ifdef verbose
  // print out the spline evaluation?
  if( Nmeas > 0 ) {
    double t = 0.0 ;
    for( t = 0.0 ; t < x[ Nmeas - 1 ] ; t+= 0.005 ) {
      fprintf( stdout , "[spline] %g %g \n" , t , 
	       cubic_eval( x , meas , der , t , Nmeas ) ) ;
    }
  }
  #endif
  // evaluate at "scale" error flag is -1
  return solve_spline( x , meas , der , scale , change_up ) ;
}

// allocate WF temporaries
int
allocate_WF( struct wflow_temps *WF ,
	     const GLU_bool memcheap ,
	     const GLU_bool adaptive )
{
  WF->lat2 = NULL ; WF->lat3 = NULL ; WF->lat4 = NULL ; WF->Z = NULL ;
  WF->lat_two = NULL ; WF->red = NULL ;
  int FLAG = GLU_SUCCESS ;
  if( adaptive == GLU_TRUE ) {
    if( ( WF->lat_two = allocate_lat( ) ) == NULL ) {
      fprintf( stderr , "[WFLOW] adaptive allocation failure \n" ) ;
      FLAG = GLU_FAILURE ;
    }
  }
  if( memcheap == GLU_TRUE ) {
    if( ( WF->lat2 = allocate_s_site( LCU , ND , NCNC ) ) == NULL ||
	( WF->Z    = allocate_s_site( LVOLUME , ND , TRUE_HERM ) ) == NULL ||
        #ifdef IMPROVED_SMEARING
	( WF->lat3 = allocate_s_site( 2*LCU , ND , NCNC ) ) == NULL ||
	( WF->lat4 = allocate_s_site( 2*LCU , ND , NCNC ) ) == NULL ||
        #else
	( WF->lat3 = allocate_s_site( LCU , ND , NCNC ) ) == NULL ||
	( WF->lat4 = allocate_s_site( LCU , ND , NCNC ) ) == NULL ||
        #endif
	( WF->red = malloc( CLINE*Latt.Nthreads*sizeof(double)) ) == NULL ) {
      fprintf( stderr , "[WFLOW] allocation failure \n" ) ;
      FLAG = GLU_FAILURE ;
    }
  } else {
    if( ( WF->lat2 = allocate_s_site( LVOLUME , ND , NCNC ) ) == NULL ||
	( WF->Z    = allocate_s_site( LVOLUME , ND , TRUE_HERM ) ) == NULL ||
	( WF->red  = malloc( CLINE*Latt.Nthreads*sizeof(double)) ) == NULL ) {
      fprintf( stderr , "[WFLOW] temporary field allocation failure\n" ) ;
      FLAG = GLU_FAILURE ;
    }
  }
  return FLAG ;
}

// free the WF temporaries
void
free_WF( struct wflow_temps *WF ,
         const GLU_bool memcheap ,
	 const GLU_bool adaptive )
{
  if( adaptive == GLU_TRUE ) {
    free_lat( WF->lat_two ) ;
  }
  if( WF -> red != NULL ) {
    free( WF->red ) ;
  }
  if( memcheap == GLU_TRUE ) {
    // free our fields
    free_s_site( WF->Z , LVOLUME , ND ) ;
#ifdef IMPROVED_SMEARING
    free_s_site( WF->lat2 , 2*LCU , ND ) ;
    free_s_site( WF->lat3 , 2*LCU , ND ) ;
    free_s_site( WF->lat4 , 2*LCU , ND ) ;
#else
    free_s_site( WF->lat2 , LCU , ND ) ;
    free_s_site( WF->lat3 , LCU , ND ) ;
    free_s_site( WF->lat4 , LCU , ND ) ;
#endif
  } else {
    // free our fields
    free_s_site( WF->lat2 , LVOLUME , ND ) ;
    free_s_site( WF->Z , LVOLUME , ND ) ;
  }
  return ;
}

// print out the general flow measurements
void
print_flow( const struct wfmeas *curr ,
	    const double err ,
	    const double delta_t)
{
  if( delta_t < 0 ) {
    fprintf( stdout , "[WFLOW-TSTOP] {err} %1.3e {t} %f {dt} %g " ,
	     err , curr -> time , delta_t ) ;
  } else {
    fprintf( stdout , "[WFLOW] {err} %1.3e {t} %f {dt} %g " ,
	     err , curr -> time , delta_t ) ;
  }
  fprintf( stdout , "{p} %g {q} %g {Gt} %g {GtP} %g\n" ,
	   curr -> avplaq , curr -> qtop , curr -> Gt ,
	   curr->time*curr->time*32*(1-curr -> avplaq) 
	   ) ;
  return ;
}

// print out the general beginning information
void
print_GG_info( void ) 
{
  fprintf( stdout , "[WFLOW] Taking ({W},{GG} and {Qtop}) measurements"
	   "from t >= %g \n" , (double)WFLOW_MEAS_START ) ;
#ifndef WFLOW_TIME_ONLY
  fprintf( stdout , "[WFLOW] fine measurements at t_0 >= %g \n" ,
	   (double)T0_STOP ) ; 
  fprintf( stdout , "[WFLOW] fine measurements at w_0 >= %g \n" ,
	   (double)W0_STOP ) ;
#endif
  fprintf( stdout , "[WFLOW] Stopping flow integration at t >= %g \n\n" , 
	   (double)WFLOW_TIME_STOP ) ; 
  return ;
}

// use the flow of G and W for scale setting
int
scaleset( struct wfmeas *curr ,
	  const size_t NT0 ,
	  const double T_0[NT0] ,
	  const size_t NW0 ,
	  const double W_0[NW0] ,
	  const size_t count ,
	  const GLU_bool is_clover ) 
{
  struct wfmeas *Head = curr ;
  
  // now we have the number of measurements in count
  double *GT   = malloc( ( count + 1 ) * sizeof( double ) ) ;
  double *time = malloc( ( count + 1 ) * sizeof( double ) ) ;
  double *der  = malloc( ( count + 1 ) * sizeof( double ) ) ;
  
  // error handling
  double t0 = -1 , w0 = -1 ;
  int flag = GLU_SUCCESS ;
  char clovstr[ 64 ] = "Plaq" ;
  if( is_clover == GLU_TRUE ) {
    sprintf( clovstr , "Clover" ) ;
  }  
  // traverse the list backwards setting the time and GT
  size_t i ;
  // traverse back down the linked list
  for( i = 0 ; i < ( count + 1 ) ; i++ ) {
    time[ count - i ] = curr -> time ;
    GT[ count - i ] = curr->time*curr->time*32*(1-curr -> avplaq) ;
    if( is_clover == GLU_TRUE ) {
      GT[ count - i ] = curr -> Gt ;
    }
    curr = curr -> next ;
  }
  if( count > 0 ) {
    for( size_t tm = 0 ; tm < NT0 ; tm++ ) {
      t0 = evaluate_scale( der , time , GT , count , T_0[tm] ) ;
      if( t0 == -1 ) {
        #pragma omp master
	{
	  fprintf( stderr , "[WFLOW] cannot compute %s t0 as we do not bracket GT\n" , clovstr ) ;
	}
	flag = GLU_FAILURE ;
	goto free ;
      }
      #pragma omp master
      {
	fprintf( stdout , "[GT-scale %s] G(%g) %1.12e \n" , clovstr , T_0[tm] , sqrt( t0 ) ) ;
      }
    }
  }
  // W(t) = t ( dG(t) / dt )
  if( count > 0 ) {
    for( i = 0 ; i < ( count + 1 ) ; i++ ) {
      GT[ i ] = time[ i ] * der[ i ] ;
    }
    for( size_t wm = 0 ; wm < NW0 ; wm++ ) {
      w0 = evaluate_scale( der , time , GT , count , W_0[wm] ) ;
      if( w0 == -1 ) {
        #pragma omp master
	{
	  fprintf( stderr , "[WFLOW] cannot compute %s w0 as we do not bracket WT\n" , clovstr ) ;
	}
	flag = GLU_FAILURE ;
	goto free ;
      }
      #pragma omp master
      {
	fprintf( stdout , "[WT-scale %s] W(%g) %1.12e \n" , clovstr , W_0[wm] , sqrt( w0 ) ) ;
      }
    }
  }
 free :
  free( der ) ;
  free( time ) ;
  free( GT ) ;

  // just reset the linked list here
  curr = Head ;

  return flag ;
}

// driver function for the more memory expensive smearing
void
step_distance( struct site *__restrict lat ,
	       struct wflow_temps WF ,
	       const double delta_t ,
	       const smearing_types SM_TYPE ,
	       void (*project)( GLU_complex log[ NCNC ] , 
				GLU_complex *__restrict staple , 
				const GLU_complex link[ NCNC ] , 
				const double smear_alpha ) )
{
  // RK3 parameters
  const double rk1 = -0.52941176470588235294*delta_t ;
  const double rk2 =  delta_t ;
  const double rk3 = -delta_t ;
  const double mthreeOfour = -3.0/4.0 ;
  const double mseventeenOthsix = -17.0/36.0 ;
  const double eightOnine = 8.0/9.0 ;
  // set z to zero
  size_t i , mu ;
#pragma omp for private(i) collapse(2)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      memset( WF.Z[i].O[mu] , 0.0 , TRUE_HERM*sizeof( GLU_complex ) ) ;
    }
  }
  // flow forwards one fictitious timestep
  RK3step( WF , lat , mseventeenOthsix , rk1 , SM_TYPE , project ) ;
  RK3step( WF , lat , eightOnine       , rk2 , SM_TYPE , project ) ;
  RK3step( WF , lat , mthreeOfour      , rk3 , SM_TYPE , project ) ;
  return ;
}

// perform a usual step
void
step_distance_memcheap( struct site *__restrict lat ,
			struct wflow_temps WF ,
			const double delta_t ,
			const smearing_types SM_TYPE ,
			void (*project)( GLU_complex log[ NCNC ] , 
					 GLU_complex *__restrict staple , 
					 const GLU_complex link[ NCNC ] , 
					 const double smear_alpha ) )
{
  // RK3 parameters
  const double rk1 = -0.52941176470588235294 * delta_t ;
  const double rk2 =  delta_t ;
  const double rk3 = -delta_t ;
  const double mthreeOfour = -3.0/4.0 ;
  const double mseventeenOthsix = -17.0/36.0 ;
  const double eightOnine = 8.0/9.0 ;
  // set z to zero
  size_t i , mu ;
#pragma omp for private(i) collapse(2)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      memset( WF.Z[i].O[mu] , 0.0 , TRUE_HERM*sizeof( GLU_complex ) ) ;
    }
  }
  // flow forwards one timestep
  RK3step_memcheap( WF , lat , mseventeenOthsix , rk1 , SM_TYPE , project ) ;
  RK3step_memcheap( WF , lat , eightOnine       , rk2 , SM_TYPE , project ) ;
  RK3step_memcheap( WF , lat , mthreeOfour      , rk3 , SM_TYPE , project ) ;
  return ;
}

// updates the measurement linked list
void
update_meas_list( struct wfmeas *head ,
		  struct wfmeas *curr ,
		  double *red ,
		  const double new_plaq ,
		  const double t ,
		  const double delta_t ,
		  const double errmax ,
		  const struct site *lat )
{
  // compute field strength tensor
  compute_Gmunu_th( red , lat ) ;
  double GT = 0.0 , Q = 0.0 ;
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    GT += red[ 0 + CLINE*i ] ;
    Q  += red[ 1 + CLINE*i ] ;
  }
  GT /= 16*LVOLUME ;
  Q *= -0.001583143494411527678811 ;
  
  curr -> time   = t ;
  curr -> Gt     = curr -> time * curr -> time * GT ;
  curr -> qtop   = Q ;
  curr -> avplaq = new_plaq ; 
#pragma omp master
  {
    print_flow( curr , errmax , delta_t ) ;
  }
  curr -> next = head ;
  return ; 
}

