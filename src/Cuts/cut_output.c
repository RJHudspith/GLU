/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (cut_output.c) is part of GLU.

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
   @file cut_output.c
   @brief output functions for the cutting procedures
 */
#include "Mainfile.h"

#include "geometry.h"  // geometry for computing momenta
#include "GLU_bswap.h" // byte swapping
#include "str_stuff.h" // append_char()

// big endian file writers, set up to mimic the c-standard fwrite but cannot
// use a multiple of the array length for the size value
static int
FWRITE( void *arr ,
	const size_t size , // MUST BE sizeof(int) or sizeof(double)
	const size_t stride ,
	FILE *file )
{
  // I enforce the re byte-swap in case we want to look at the data again
  if( size == sizeof( int ) ) {
    if( !WORDS_BIGENDIAN ) { bswap_32( stride , arr ) ; }
    fwrite( arr , size , stride , file ) ;
    if( !WORDS_BIGENDIAN ) { bswap_32( stride , arr ) ; }
  } else if ( size == sizeof( double ) ) {
    if( !WORDS_BIGENDIAN ) { bswap_64( stride , arr ) ; }
    fwrite( arr , size , stride , file ) ;
    if( !WORDS_BIGENDIAN ) { bswap_64( stride , arr ) ; }
  } else {
    fprintf( stderr , "[CUT] I do not understand your byte-swapping "
	     "size %zu \n" , size ) ;
    return GLU_FAILURE ;
  }
  return GLU_SUCCESS ;
}

// gives us some nominal information on the cutting procedure
int
check_psq( const struct cut_info CUTINFO )
{
  size_t small = Latt.dims[0] ;
  size_t i ;
  for( i = 1 ; i < ND ; i++ ) {
    if( Latt.dims[i] < small ) {
      small = Latt.dims[i] ;
    }
  }

  // Type is the type of cut we are using
  // 0 is normal hypercubic
  // 1 is normal psq
  // 2 is a cylinder cut
  // 3 is with an angle

  fprintf( stdout , "\n[CUTS] " ) ;
  if( CUTINFO.dir == NONEXCEPTIONAL ) {
    fprintf( stdout , "PSQ Cutting Procedure \n" ) ;
  } else {
    switch( CUTINFO.type ) {
    case HYPERCUBIC_CUT :
      fprintf( stdout , "HYPERCUBIC Cutting Procedure \n" ) ;
      break ;
    case PSQ_CUT :
      fprintf( stdout , "PSQ Cutting Procedure \n" ) ;
      break ;
    case CYLINDER_CUT :
      fprintf( stdout , "Cylinder Cutting Procedure \n" ) ;
      fprintf( stdout , "[CUTS] Cylinder Width in Lattice units :: %g \n" , 
	       CUTINFO.cyl_width ) ; 
      fprintf( stdout , "[CUTS] Physical Momentum Cap :: %g \n" , 
	       TWOPI*CUTINFO.cyl_width/small  ) ; 
      break ;
    case CYLINDER_AND_CONICAL_CUT :
      fprintf( stdout , "Cylinder and Conical Cutting Procedure \n" ) ;
      fprintf( stdout , "[CUTS] Cylinder Radius in Lattice units :: %g \n" , 
	       CUTINFO.cyl_width ) ; 
      fprintf( stdout , "[CUTS] Physical Momentum Cap :: %g \n" , 
	       TWOPI*CUTINFO.cyl_width/small ) ; 
      fprintf( stdout , "[CUTS] Cone Angle :: %g \n" , 
	       (double)CUTINFO.angle ) ; 
      break ;
    default:
      fprintf( stderr , "I don't recognise the type.. Leaving \n" ) ;
      return GLU_FAILURE ;
    }
  }
  fprintf( stdout , "[CUTS] Maximum p^2 :: %zu \n\n" , CUTINFO.max_mom ) ;
  if( CUTINFO.definition == LOG_DEF ) {
    fprintf( stdout , "[CUTS] Using the LOGARITHMIC field definition \n" ) ;
  } else {
    fprintf( stdout , "[CUTS] Using the LINEAR field definition \n" ) ;
  }
  return GLU_SUCCESS ;
}

// Decides the shape of our output string
char *
output_str_struct( const struct cut_info CUTINFO )
{
  char tmp[ 256 ] ;

  #ifndef CONDOR_MODE
  char *str = malloc( (1+strlen(CUTINFO.where))*sizeof( char ) ) ;
  sprintf( str , "%s" , CUTINFO.where ) ;
  fprintf( stdout , "[CUTS] Our file is located @ %s \n", str ) ;
  #else
  char *str = malloc( (1+strlen("./"))*sizeof( char ) ) ;
  sprintf( str , "./" ) ;
  #endif

  // each mode has a slightly different start so that you know what you are getting into 
  switch( CUTINFO.dir ) {
  case INSTANTANEOUS_GLUONS :
    append_char( &str , "G2spG2t" ) ;
    break ;
  case SMEARED_GLUONS :
    append_char( &str , "G2_G2SM" ) ;
    break ;
  case GLUON_PROPS :
    append_char( &str , "G2F2" ) ;
    break ;
  case CONFIGSPACE_GLUONS :
    append_char( &str , "CSPACE" ) ;
    break ;
  case FIELDS :
    append_char( &str , "MATS" ) ;
    break ;
  case EXCEPTIONAL :
    append_char( &str , "MOMgg" ) ;
    break ;
  case NONEXCEPTIONAL :
#ifdef PROJ_GRACEY
    append_char( &str , "spgracey" ) ;
#else
    append_char( &str , "spboucaud" ) ;
#endif
    break ;
  case STATIC_POTENTIAL :
    append_char( &str , "STATPOT" ) ;
    break ;
  case TOPOLOGICAL_SUSCEPTIBILITY :
    append_char( &str , "QSUSC" ) ;
    break ;
  default : // default behaviour is the 3 gluon vertex
    append_char( &str , "MOMgg" ) ;
    break ;
  }
  // color info
  sprintf( tmp , "_SU%d_" , NC ) ;
  append_char( &str , tmp ) ;

  // print out the dimensions of the lattice used for the cut
  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    sprintf( tmp , "%zux" , Latt.dims[mu] ) ;
    append_char( &str , tmp ) ;
  }
  sprintf( tmp , "%zu_" , Latt.dims[mu] ) ;
  append_char( &str , tmp ) ;

  if( CUTINFO.dir == CONFIGSPACE_GLUONS ) {
    if( CUTINFO.definition == LOG_DEF ) {
      sprintf( tmp , "LOG.%zu.bin" , Latt.flow ) ;
      append_char( &str , tmp ) ;
    } else {
      sprintf( tmp , "LIN.%zu.bin" , Latt.flow ) ;
      append_char( &str , tmp ) ;
    }
    return str ;
  }

  /////// MOMENTUM SPACE STUFF FROM HERE ///////
  sprintf( tmp , "%zu_" , CUTINFO.max_mom ) ;
  append_char( &str , tmp ) ;

  // static potential measurement has T-detail here
  if( CUTINFO.dir == STATIC_POTENTIAL ) { 
    sprintf( tmp , "T%zu_" , CUTINFO.max_t ) ;
    append_char( &str , tmp ) ;
  }

  // switch on the correct cut, only PSQ is allowed for the non-exceptional
  momentum_cut_def type = CUTINFO.type ;
  if( CUTINFO.dir == NONEXCEPTIONAL ) {
    type = PSQ_CUT ;
  }

  // print out the cutting type performed
  switch( type ) {
  case HYPERCUBIC_CUT :
    append_char( &str , "cutHYP" ) ;
    break;
  case PSQ_CUT :
    append_char( &str , "cutPSQ" ) ;
    break;
  case CYLINDER_CUT :
    sprintf( tmp , "cutCYL_w%g" , CUTINFO.cyl_width ) ;
    append_char( &str , tmp ) ;
    break ;
  case CYLINDER_AND_CONICAL_CUT :
    sprintf( tmp , "cutCON_%zu" , CUTINFO.angle ) ;
    append_char( &str , tmp ) ;
    break ; 
  default:
    break ; 
  }

  // print out the momentum def used for the projectors, only used for
  // the three point functions
  if( CUTINFO.dir != STATIC_POTENTIAL && 
      CUTINFO.dir != TOPOLOGICAL_SUSCEPTIBILITY ) {
    #ifdef SIN_MOM
    append_char( &str , "_SINMOM" ) ;
    #else
    append_char( &str , "_PSQMOM" ) ;
    #endif
  }

  // have PROJ_GRACEY == 0 for the 6 combined
  #ifdef PROJ_GRACEY
  if( CUTINFO.dir == NONEXCEPTIONAL ) {
    #if PROJ_GRACEY > -1 && PROJ_GRACEY < 15
    sprintf( tmp , "_PROJ%d" , PROJ_GRACEY ) ;
    append_char( &str , tmp ) ;
    #endif
  }
  #endif
  
  // field definition for momentum space glue
  if( CUTINFO.dir != STATIC_POTENTIAL && 
      CUTINFO.dir != TOPOLOGICAL_SUSCEPTIBILITY ) {
    if( CUTINFO.definition == LOG_DEF ) {
      append_char( &str , "_LOG" ) ;
    } else {
      append_char( &str , "_LIN" ) ;
    }
  }
  // and finish up with the configuration number
  sprintf( tmp , ".%zu.bin" , Latt.flow ) ;

  append_char( &str , tmp ) ;

  return str ;
}

// Write the gluonic two and three point functions to the file Ap
void
write_complex_g2g3_to_list( FILE *Ap , 
			    double complex *g2 , 
			    double complex *g3 , 
			    size_t num_mom[ 1 ] ) 
{
#ifdef ASCII_CHECK
  printf("%d \n", num_mom[ 0 ] ) ;
  size_t i ;
  for( i = 0 ; i < num_mom[0] ; i++ ) { printf("%e \n", g2[i] ) ; }
  printf("%d \n", num_mom[ 0 ] ) ;
  for( i = 0 ; i < num_mom[0] ; i++ ) { printf("%e \n", g3[i] ) ; }
#else
  uint32_t nmom[1] = { num_mom[0] } ;
  // write NUM_MOM[0] , G2 , NUM_MOM[0] , G3
  FWRITE( nmom , sizeof(uint32_t) , 1 , Ap ) ;
  FWRITE( g2 , sizeof(double) , 2 * num_mom[0] , Ap ) ;
  FWRITE( nmom , sizeof(uint32_t) , 1 , Ap ) ;
  FWRITE( g3 , sizeof(double) , 2 * num_mom[0] , Ap ) ;
#endif
  return ;
}

// Write the gluonic two point function to the file Ap
void
write_g2_to_list( FILE *Ap , 
		  double *g2 , 
		  size_t num_mom[ 1 ] ) 
{
#ifdef ASCII_CHECK
  size_t i ;
  fprintf( stdout , "[CUTS] G2\n%d\n" , num_mom[0] ) ;
  for( i = 0 ; i < num_mom[0] ; i++ ) { 
    fprintf( stdout , "%e\n", g2[i] ) ; 
  }
  fprintf( stdout , "%d\n", num_mom[ 0 ] ) ;
#else
  // write NUM_MOM[0] , G2
  uint32_t nmom[ 1 ] = { num_mom[0] } ;
  FWRITE( nmom , sizeof(uint32_t) , 1 , Ap ) ;
  FWRITE( g2 , sizeof(double) , num_mom[0] , Ap ) ;
#endif
  return ;
}

// Write the gluonic two and three point functions to the file Ap
void
write_g2g3_to_list( FILE *Ap , 
		    double *g2 , 
		    double *g3 , 
		    size_t num_mom[ 1 ] ) 
{
#ifdef ASCII_CHECK
  fprintf( stout , "[CUTS] g2g3\n%d\n", num_mom[ 0 ] ) ;
  size_t i ;
  for( i = 0 ; i < num_mom[0] ; i++ ) { 
    printf( stdout , "%e\n", g2[i] ) ; 
  }
  printf("%d \n", num_mom[ 0 ] ) ;
  for( i = 0 ; i < num_mom[0] ; i++ ) { 
    printf( stdout , "%e\n", g3[i] ) ; 
  }
#else
  uint32_t nmom[1] = { num_mom[0] } ;
  // write NUM_MOM[0] , G2 , NUM_MOM[0] , G3
  FWRITE( nmom , sizeof(uint32_t) , 1 , Ap ) ;
  FWRITE( g2 , sizeof(double) , nmom[0] , Ap ) ;
  FWRITE( nmom , sizeof(uint32_t) , 1 , Ap ) ;
  FWRITE( g3 , sizeof(double) , nmom[0] , Ap ) ;
#endif
  return ;
}

void
write_lattice_fields( FILE *Ap , 
		      const struct site *A , 
		      const struct veclist *list , 
		      size_t num_mom[1] )
{
  uint32_t nmom[1] = { num_mom[0] } ;
  // write the length of the array
  FWRITE( nmom , sizeof(uint32_t) , 1 , Ap ) ;

  //write lattice fields
  const size_t stride = ND * NCNC * 2 ; // the 2 is for the re/im parts
  double *aout = malloc( stride * sizeof( double ) ) ;

  size_t i ;
  for( i = 0 ; i < num_mom[0] ; i++ ) {
    size_t mu ;
    #pragma omp parallel for private(mu)
    for( mu = 0 ; mu < ND ; mu++ ) {
      size_t a , b ;
      for( a = 0 ; a < NCNC ; a++ ) {
	b = 2 * ( a + mu * NCNC ) ;
	aout[ b ] = (double)creal( A[ list[i].idx ].O[mu][a] ) ;
	aout[ b + 1 ] = (double)cimag( A[ list[i].idx ].O[mu][a] ) ;
      }
    }
    FWRITE( aout , sizeof(double) , stride , Ap ) ;
  }
  free( aout ) ;
  return ;
}

// Support for writing our momentum list
void
write_mom_veclist( FILE *Ap , 
		   size_t *num_mom , 
		   const struct veclist *list ,
		   const size_t DIR )
{
  const int stride = ( DIR + 1 ) * num_mom[ 0 ] ;
  int *kall = malloc( stride * sizeof(int) ) ;

  uint32_t nmom[1] = { num_mom[0] } ;
  FWRITE( nmom , sizeof(int) , 1 , Ap ) ;

  size_t i ;
  //write momenta
#pragma omp parallel for private(i)
  for( i = 0 ; i < num_mom[0] ; i++ ) {
    size_t n = i * ( DIR + 1 ) ;
    kall[ n ] = DIR ; 
    n++ ; 
    size_t a ;
    for( a = 0  ;  a < DIR  ;  a++ ) {
      kall[ n ] = list[i].MOM[a] ; 
      n++ ; 
    }
  }
  FWRITE( kall , sizeof(int) , stride , Ap ) ;
  free( kall ) ;
  return ;
}

// write out some moments, re and imaginary
void
write_moments( struct Qmoments *Qmom ,
	       const size_t Nmeasurements )
{
  FILE *Qm = NULL , *Q2m = NULL ;
  double *re = malloc( NQMOMENTS * sizeof( double ) ) ;
  uint32_t nmom[ 1 ] = { NQMOMENTS } ;
  uint32_t Nmeas[ 1 ] = { Nmeasurements } ;
  char filename[ 256 ] ;
  
  size_t mu , i , j ;
  for( mu = 0 ; mu < ND ; mu++ ) {

    sprintf( filename , "Qmoment_d%zu.%zu.bin" , mu , Latt.flow ) ;
    Qm = fopen( filename , "wb" ) ;
    sprintf( filename , "Q2moment_d%zu.%zu.bin" , mu , Latt.flow ) ;
    Q2m = fopen( filename , "wb" ) ;

    FWRITE( Nmeas , sizeof( uint32_t ) , 1 , Qm ) ;
    FWRITE( Nmeas , sizeof( uint32_t ) , 1 , Q2m ) ;
    
    for( i = 0 ; i < Nmeasurements ; i++ ) {
      for( j = 0 ; j < NQMOMENTS ; j++ ) {
	re[j] = Qmom[i].Q[ j+mu*NQMOMENTS ] ;
      }
      FWRITE( nmom , sizeof( uint32_t ) , 1 , Qm ) ;
      FWRITE( re , sizeof( double ) , NQMOMENTS , Qm ) ;
      
      for( j = 0 ; j < NQMOMENTS ; j++ ) {
	re[j] = Qmom[i].Q2[j + mu*NQMOMENTS] ;
      }
      FWRITE( nmom , sizeof( uint32_t ) , 1 , Q2m ) ;
      FWRITE( re , sizeof( double ) , NQMOMENTS , Q2m ) ;
    }

    fclose( Qm ) ; fclose( Q2m ) ; 
  }
  
  free( re ) ;
  
  return ;
}

// write out rsq's for the correlator
void
write_rr_values( FILE *Ap ,
		 size_t size[1] ,
		 const size_t *rsq ,
		 const size_t max_r2 ,
		 const size_t ARR_SIZE )
{
  size_t i ;
#ifdef ASCII_CHECK
  fprintf( stdout , "[CUTS] rr\n%d\n" , size[0] ) ;
  // loop the possible rsq's
  for( i = 0 ; i < ARR_SIZE ; i++ ) {
    if( rsq[i] < max_r2 ) {
      size_t rr[1] = { rsq[i] } ;
      printf( "%zu\n" , rr[0] ) ;
    }
  }
#endif
  uint32_t nmom[1] = { size[0] } ;
  // output to a file ... for reading in with UKhadron or something
  // write out the length of the array ...
  FWRITE( nmom , sizeof( uint32_t ) , 1 , Ap ) ;
  // loop the possible rsq's
  for( i = 0 ; i < ARR_SIZE ; i++ ) {
    if( rsq[i] < max_r2 ) {
      uint32_t rr[1] = { rsq[i] } ;
      FWRITE( rr , sizeof(uint32_t) , 1 , Ap ) ;
    }
  }
}

// Support for writing our triplets momentum
void
write_triplet_mom_list( FILE *Ap , 
			size_t *num_mom , 
			int **momentum ,
			int **triplet )
{
  const size_t stride = ( ND + 1 ) * num_mom[0] ;
  int *kall = malloc( stride * sizeof( int ) ) ;

#ifdef ASCII_CHECK
  fprintf( stdout , "[CUTS] triplet momlist\n%d\n", num_mom[0] ) ;
#endif

  uint32_t nmom[1] = { num_mom[0] } ;
  FWRITE( nmom , sizeof(int) , 1 , Ap ) ;

  // start writing the momenta to the list, we just write the first momentum of the triplet...
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < num_mom[0] ; i++ ) {
    size_t n = i * ( ND + 1 ) , mu ;
    kall[ n ] = ND ; 
    n++ ; 
    for( mu = 0 ; mu < ND ; mu++ ) {
      kall[n] = momentum[ triplet[ i ][ 0 ] ][ mu ] ;
      n++ ;
    }
  }

#ifdef ASCII_CHECK
  size_t mu , b = 0 ;
  for( i = 0 ; i < num_mom[0] ; i++ ) {
    for( mu = 0 ; mu < ND + 1 ; mu++ ) {
      fprintf( stdout , "%d ", kall[ b ] ) ; 
      b++ ;
    }
    fprintf( stdout , "\n" ) ;
  }
  fprintf( stdout , "%d \n", kall[b] ) ;
  fprintf( stdout , "%zu \n", num_mom[0] ) ;
#endif
  FWRITE( kall , sizeof(int) , stride , Ap ) ;
  free( kall ) ;
  return ;
}

// write the timeslices, LT[0] is the size of the array
void
write_tslice_list( FILE *Ap , 
		   size_t *LT )
{
  size_t t ;
  uint32_t T[ LT[0] ] , lt[1] = { *LT } ;
#ifdef ASCII_CHECK
  fprintf( stdout , "[CUTS] tslice list\n%d\n" , LT[0] ) ;
  for( t = 0 ; t < LT[0] ; t++ ) { 
    fprintf( stdout , "%d\n" , t ) ; 
  }
  fprintf( stdout , "%d\n" , LT[0] ) ;
#endif
  FWRITE( lt , sizeof(uint32_t) , 1 , Ap ) ;
  for( t = 0 ; t < LT[0] ; t++ ) { T[ t ] = t ; }
  FWRITE( T , sizeof(uint32_t) , LT[0] , Ap ) ;
  return ;
}

