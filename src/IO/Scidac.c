/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (Scidac.c) is part of GLU.

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
   @file Scidac.c
   @brief Scidac and ILDG configuration header readers and writers
 */
#include "Mainfile.h"

#include "GLU_bswap.h" // for the byteswaps
#include "GLU_timer.h" // for the date
#include "str_stuff.h" // append_char()
#include "XML_info.h"  // gets the important info from the xml header info

// < Message Begin Mask (Internal) 
#define MB_MASK ((unsigned char)0x80)

// < Message End Mask (Internal)
#define ME_MASK ((unsigned char)0x40)

// maximum header elements
#define MAX_HDR64 18

// end of file flag
#define GLU_EOF 1 // GLU_SUCCESS is 0, GLU_FAILURE is -1

// this contains all the header data
static union {
  uint64_t int64[ MAX_HDR64 ] ;
  uint32_t int32[ 2*MAX_HDR64 ] ;
  uint16_t int16[ 4*MAX_HDR64 ] ;
  unsigned char uchr[ 8*MAX_HDR64 ] ;
} header ;

// some padding
static inline int 
lime_padding( const size_t nbytes ){ return ( nbytes%8 != 0 ) ? 8 - (nbytes%8) : 0 ; }
// accessors and what have you for the header union
static inline uint64_t
header_datalength( void ) { return header.int64[1] ; }
static inline uint32_t
magic_number( void ) { return header.int32[0] ; }
static inline uint16_t
header_version( void ) { return header.int16[2] ; }

static unsigned char *lime_hdr_rec_type = &header.uchr[16] ;

/**
   @fn static int parse_SCIDAC_hdr( FILE *infile , struct head_data *HEAD_DATA , const int record )
   @brief looks through the scidac header setting up file IO stuff
 */
static int
parse_SCIDAC_hdr( FILE *infile , 
		  struct head_data *HEAD_DATA )
{
  // assumes we have opened it fine and are at the start of the file
  const char myname[] = "GLU::parse_SCIDAC_hdr";

  // fread the whole header ...
  const int status = fread( (void*)header.int64 , sizeof( int64_t ) , 
			    MAX_HDR64 , infile ) ;
  if( status != MAX_HDR64 ) { return GLU_EOF ; }
  
  uint32_t i_magic_no[1] = { magic_number( ) } ;
  if( !WORDS_BIGENDIAN ) { bswap_32( 1 , i_magic_no ) ; }

  if( i_magic_no[0] != 1164413355 ) {
    fprintf( stderr , "%s wrong magic number %x vs. %x \n" , myname ,
	    i_magic_no[0] , 1164413355 ) ;
    return GLU_FAILURE ;
  }

  // now we look at the version number
  unsigned int i_version[1] = { header_version( ) } ;
  if( !WORDS_BIGENDIAN ) { bswap_32( 1 , i_version ) ; }

  uint64_t i_data_length[1] = { header_datalength() } ;
  if( !WORDS_BIGENDIAN ) { bswap_64( 1 , i_data_length ) ; }

  const size_t padding = lime_padding( (size_t)i_data_length[0] ) ;

  // again, typebuf is uninteresting ...
  unsigned char *typebuf = (unsigned char*)lime_hdr_rec_type ;
  // if we hit the binary data we exit for a binary read, could pass the length
  // of the data but the header has already given us the logical dimensions
  if( strcmp( (const char*)typebuf , "scidac-binary-data" ) == 0 ||
      strcmp( (const char*)typebuf , "ildg-binary-data" ) == 0 ) {
    // and so we exit ...
    return GLU_EOF ; // is not end of file, but not error either
  }
  // we should be able to read the remaining data
  char *io_data = malloc( ( i_data_length[0] + 1 ) * sizeof( char ) ) ;
  if( fread( io_data , sizeof( char ) , i_data_length[0] , infile ) != 
      i_data_length[0] ) {
    fprintf( stderr , "[IO] xml information reading failed ... Leaving \n" ) ;
    free( io_data ) ;
    return GLU_FAILURE ;
  }
  io_data[ i_data_length[0] ] = '\0' ;
  // BQCD put their plaquette and cksum at the end without xml tags hmpf!
  if( strcmp( (const char*)typebuf , "bqcd-plaquette" ) == 0 ) {
    struct head_data tmp = *HEAD_DATA ;
    char *endptr ;
    tmp.plaquette = strtod( io_data , &endptr ) ;
    if( io_data == endptr ) {
      printf( "[IO] BQCD plaquette misread \n" ) ;
      return GLU_FAILURE ;
    }
    *HEAD_DATA = tmp ;
  }
  // read in the possible checksum
  if( strcmp( (const char*)typebuf , "bqcd-cksum" ) == 0 ) {
    struct head_data tmp = *HEAD_DATA ;
    unsigned long int cksum ;
    sscanf( io_data , "%lu" , &cksum ) ;
    tmp.checksum = (uint32_t)cksum ;
    *HEAD_DATA = tmp ;
  }
  parse_and_set_xml_SCIDAC( io_data , HEAD_DATA ) ;

  // free for now
  free( io_data ) ;

  // and skip the file along by "padding" amount
  if( padding > 0 ) {
    char pad_str[ padding ] ;
    if( fread( pad_str , sizeof( char ) , padding , infile ) != padding ) {
      fprintf( stderr , "[IO] Reading of padded data failed \n" ) ;
      return GLU_FAILURE ;
    }
  }
  return GLU_SUCCESS ;
}

// for writing the header record
static int
write_SCIDAC_binary( FILE *__restrict out ,
		     char *record ,
		     const size_t record_length ,
		     const char *banf )
{
  // padding ...
  const unsigned char padbuf[7] = {0x00,0x00,0x00,0x00,0x00,0x00,0x00};
  // we zero the header
  int i ;
  for( i = 0 ; i < MAX_HDR64 ; i++ ) { header.int64[i] = 0 ; }
  // is magic
  uint32_t i_magic_no[1] = { 1164413355 } ;
  if( !WORDS_BIGENDIAN ) { bswap_32( 1 , i_magic_no ) ; }
  header.int32[0] = i_magic_no[0] ;
  // is version
  uint16_t i_version[1] = { 1 } ;
  if( !WORDS_BIGENDIAN ) { bswap_16( 1 , i_version ) ; }
  header.int16[2] = i_version[0] ;
  // is ME and MB
  header.uchr[6] = MB_MASK ;
  // is datalength
  uint64_t datalength[1] = { record_length } ;
  if( !WORDS_BIGENDIAN ) { bswap_64( 1 , datalength ) ; }
  header.int64[1] = datalength[0] ;
  // poke in some stuff for the header
  size_t chr , ref = 16 ; // this is where the useful info stops
  for( chr = 0 ; chr < strlen( banf ) ; chr++ ) {
    header.uchr[ ref + chr ] = banf[ chr ] ;
  }
  // and write out the record
  if( fwrite( (void*)header.int64 , sizeof( uint64_t ) , MAX_HDR64 , out ) 
      != MAX_HDR64 ) {
    fprintf( stderr , "[IO] SCIDAC header writing failure ... Leaving\n" ) ;
    return GLU_FAILURE ;
  }
  fprintf( out , "%s" , record ) ;
  // write out some padding ...
  const int padding = lime_padding( record_length ) ;
  fwrite( padbuf , sizeof( unsigned char ) , padding , out ) ;
  return GLU_SUCCESS ;
}

// scidac file has a header, plus some xml data
// all in big endian, which is nice, and is called for the ILDG too
int
get_header_data_SCIDAC( FILE *infile ,
			struct head_data *HEAD_DATA )
{
  // loop the records, zealously set max number to be 12
  int record = 0 , validity ;
  while( record < 12 ) { // hmmm
    validity = parse_SCIDAC_hdr( infile , HEAD_DATA ) ; record++ ;
    if( validity == GLU_FAILURE ) { // if the reader doesn't work
      return GLU_FAILURE ;
    } else if( validity == GLU_EOF ) {
      break ;
    }
  }
  return GLU_SUCCESS ;
}

// for writing out SCIDAC configurations ...
void
write_header_SCIDAC( FILE *__restrict out )
{
  // initialisation 
  const char *start = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" ;
  char *str = malloc( (1+strlen(start))*sizeof(char) ) ;
  sprintf( str , "%s" , start ) ;
  
  // sprintf in some xml stuff
  char tmp[128] ;
  sprintf( tmp , "<scidacFile><version>1.1</version>"
	   "<spacetime>%d</spacetime><dims>" , ND ) ;
  append_char( &str , tmp ) ;
  
  int mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    sprintf( tmp , "%zu " , Latt.dims[mu] ) ;
    append_char( &str , tmp ) ;
  }
  sprintf( tmp , "%zu" , Latt.dims[mu] ) ;
  append_char( &str , tmp ) ;
  append_char( &str , "</dims><volfmt>0</volfmt></scidacFile>" ) ;

  // some dummy file
  write_SCIDAC_binary( out , str , strlen( str ) , "scidac-private-file-xml" ) ;
  free( str ) ;
  write_SCIDAC_binary( out , "<title>GLU SCIDAC archival gauge configuration</title>" ,
		       strlen( "<title>GLU SCIDAC archival gauge configuration</title>" ) ,
		       "scidac-file-xml" ) ;
  
  // gauge file informations
  char *gstr = malloc( (1+strlen(start))*sizeof(char) ) ;
  sprintf( gstr , "%s" , start ) ;
  sprintf( tmp , "<scidacRecord><version>1.0</version><date>%s GMT</date>"
	   "<globaldata>0</globaldata>" , get_date( ) ) ;
  append_char( &gstr , tmp ) ;
#ifdef SINGLE_PREC
  sprintf( tmp , "<datatype>QDP_F%d_ColorMatrix</datatype>"
	   "<precision>F</precision>" , NC ) ;
  append_char( &gstr , tmp ) ;
#else
  sprintf( tmp , "<datatype>QDP_D%d_ColorMatrix</datatype>"
	   "<precision>D</precision>" , NC ) ;
  append_char( &gstr , tmp ) ;
#endif
  sprintf( tmp , "<colors>%d</colors><typesize>%d</typesize>"
	   "<datacount>%d</datacount></scidacRecord>" , 
	   NC ,  ND * NCNC * 2 * (int)( sizeof( GLU_real )
					/ sizeof( float ) ) , ND ) ;
  append_char( &gstr , tmp ) ;
  write_SCIDAC_binary( out , gstr , strlen( gstr ) , 
		       "scidac-private-record-xml" ) ;
  free( gstr ) ;

  // dummy again
  write_SCIDAC_binary( out , "<info>GLU library configuration file</info>" ,
		       strlen( "<info>GLU library configuration file</info>" ) ,
		       "scidac-record-xml" ) ;
  
  // and then leave space for the gauge field
  char *nada = "" ; 
  write_SCIDAC_binary( out , nada , sizeof( GLU_real ) * \
		       LVOLUME * ND * NCNC * 2 ,// the 2 is for complex
		       "scidac-binary-data" ) ; 
  return ;
}

// writes the final header for the scidac record
void
write_trailing_header_SCIDAC( FILE *__restrict out ,
			      const uint32_t cksuma ,
			      const uint32_t cksumb )
{
  const char *start = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><scidacChecksum>" ;
  char *str = malloc( (1+strlen(start))*sizeof(char) ) , tmp[128] ;
  sprintf( str , "%s" , start ) ; 
  sprintf( tmp , "<version>1.0</version><suma>%x</suma><sumb>%x</sumb>" ,
	   cksuma , cksumb ) ;
  append_char( &str , tmp ) ;
  append_char( &str , "</scidacChecksum>" ) ;
  write_SCIDAC_binary( out , str , strlen( str ) , "scidac-checksum" ) ;
  free( str ) ;
  return ;
}

// and a header writer for the ILDG configuration
// for writing out SCIDAC configurations ...
void
write_header_ILDG( FILE *__restrict out )
{
  // initialisation 
  const char *start = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" ;
  char *str = malloc( (1+strlen(start))*sizeof( char ) ) , tmp[ 128 ] ;
  sprintf( str , "%s" , start ) ;
  // sprintf in some xml stuff
  sprintf( tmp , "<scidacFile><version>1.1</version>"
	   "<spacetime>%d</spacetime><dims>" , ND ) ;
  append_char( &str , tmp ) ;
  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    sprintf( tmp , "%zu " , Latt.dims[mu] ) ;
    append_char( &str , tmp ) ;
  }
  sprintf( tmp , "%zu" , Latt.dims[mu] ) ;
  append_char( &str , tmp ) ;

  append_char( &str , "</dims><volfmt>0</volfmt></ScidacFile>" ) ;
  // some dummy file
  write_SCIDAC_binary( out , str , strlen( str ) , "scidac-private-file-xml" ) ;
  free( str ) ;
  write_SCIDAC_binary( out , "<title>GLU ILDG archival gauge configuration</title>" ,
		       strlen( "<title>GLU ILDG archival gauge configuration</title>" ) ,
		       "scidac-file-xml" ) ;

  // SCIDAC gauge file informations
  char *gstr = malloc( (1+strlen(start))*sizeof( char ) ) ;
  sprintf( gstr , "%s" , start ) ;
  sprintf( tmp , "<scidacRecord><version>1.0</version><date>%s GMT</date>"
	   "<globaldata>0</globaldata>" , get_date( ) ) ;
  append_char( &gstr , tmp ) ;
#ifdef SINGLE_PREC
  sprintf( tmp , "<datatype>QDP_F%d_ColorMatrix</datatype>"
	   "<precision>F</precision>" , NC ) ;
  append_char( &gstr , tmp ) ;
#else
  sprintf( tmp , "<datatype>QDP_D%d_ColorMatrix</datatype>"
	   "<precision>D</precision>" , NC ) ;
  append_char( &gstr , tmp ) ;
#endif
  sprintf( tmp , "<colors>%d</colors><typesize>%d</typesize>"
	   "<datacount>%d</datacount></scidacRecord>" , 
	   NC , ND * NCNC * 2 * (int)( sizeof( GLU_real )
				       / sizeof( float ) ) , ND ) ;
  append_char( &gstr , tmp ) ;
  write_SCIDAC_binary( out , gstr , strlen( gstr ) ,
		       "scidac-private-record-xml" ) ;
  free( gstr ) ;
  // and finally an xml
  write_SCIDAC_binary( out , "<info>GLU library configuration file</info>" ,
		       strlen( "<info>GLU library configuration file</info>" ) ,
		       "scidac-record-xml" ) ;

  // FINALLY WE HAVE THE ILDG FORMAT STUFF //
  const char *ILDG_start = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
    "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\" xmlns:xsi=\""
    "http://www.w3.org/2001/XMLSchema-instance\" xsi:"
    "schemaLocation=\"http://www.lqcd.org/ildg/filefmt.xsd\">"
    "<version>1.0</version>" ;
  const char *ILDG_end = "</ildgFormat>" ;
  char *ILDG_str = malloc( (1+strlen(ILDG_start))*sizeof(char) ) ;
  // begin with the usual spiel and the versioning ...
  sprintf( ILDG_str , "%s" , ILDG_start ) ;
#ifdef SINGLE_PREC
  append_char( &ILDG_str , "<field>su3gauge</field><precision>32</precision>" ) ;
#else
  append_char( &ILDG_str , "<field>su3gauge</field><precision>64</precision>" ) ;
#endif
  // and input the geometry .. I enumerate the extra (>4) dimensions, always leaving t running slowest
  const char open[4][5] = { "<lx>" , "<ly>" , "<lz>" , "<lt>" } ;
  const char close[4][6] = { "</lx>" , "</ly>" , "</lz>" , "</lt>" } ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    if( mu < 4 ) {
      sprintf( tmp , "%s%zu%s" , open[mu] , Latt.dims[mu] , close[mu] ) ;
      append_char( &ILDG_str , tmp ) ;
    } else {
      sprintf( tmp , "<l%zu>%zu</l%zu>" , mu , Latt.dims[mu] , mu ) ;
      append_char( &ILDG_str , tmp ) ;
    }
  }
  sprintf( tmp , "%s%zu%s%s" , open[3] , Latt.dims[ND-1] , close[3] , ILDG_end ) ;
  append_char( &ILDG_str , tmp ) ;
  write_SCIDAC_binary( out , ILDG_str , strlen( ILDG_str ) , "ildg-format" ) ;
  free( ILDG_str ) ;

  // Logical file name; just leave it blank for now
  write_SCIDAC_binary( out , "lfn://" , strlen( "lfn://" ) , "ildg-data-lfn" ) ;
	   
  // and then leave space for the gauge field
  char *nada = "" ; 
  write_SCIDAC_binary( out , nada , sizeof( GLU_real )*LVOLUME*ND*NCNC*2 , 
		       "ildg-binary-data" ) ; // the 2 is for complex
  return ;
}

// clean this up
#ifdef DEBUG
  #undef DEBUG
#endif

// and undefs to clear up local macros
#undef A_BIG_NUMBER 
#undef MB_MASK
#undef ME_MASK
#undef MAX_HDR64
#undef GLU_EOF
