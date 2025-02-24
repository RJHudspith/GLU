/*
Copyright 2013-2025 Renwick James Hudspith

    This file (random_config.h) is part of GLU.

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
   @file random_config.h
   @brief function prototypes for lattice reunitarisation and randomisation
 */
#ifndef GLU_RANDOM_CONFIG_H
#define GLU_RANDOM_CONFIG_H

/**
   @fn void latt_reunitg( GLU_complex *__restrict *__restrict gauge )
   @brief reunitarises the gauge field
   @param gauge :: lattice gauge field
   sometimes the gauge transfomations need to be forced into SU(NC) again
 */
void 
latt_reunitg( GLU_complex *__restrict *__restrict gauge ) ;

/**
   @fn void latt_reunitU( struct site *__restrict lat )
   @brief reunitarises the gauge field
   @param lat :: lattice gauge field
  sometimes the gauge field need to be forced into SU(NC) again. Specifically for the axial gauge where the gauge transforms are nested products of what was there before can cause round-off errors
 */
void 
latt_reunitU( struct site *__restrict lat ) ;

/**
   @fn void random_gtrans( struct site *__restrict lat )
   @brief generates a random gauge transform on the field
   @param lat :: lattice gauge field
 */
void 
random_gtrans( struct site *__restrict lat ) ;

/**
   @fn void random_gtrans_slice( GLU_complex *__restrict *__restrict slice_gauge ) 
   @brief create random gauge transformation on a time-slice
   @param slice_gauge :: gauge transform matrices on a time slice
 */
void
random_gtrans_slice( GLU_complex *__restrict *__restrict slice_gauge ) ;

/**
   @fn void random_transform( struct site *__restrict lat , GLU_complex *__restrict *__restrict gauge )
   @brief compute a random gauge transform on the links.
   @param lat :: lattice gauge fields 
   @param gauge :: gauge transformation matrices
   the rationale for having this instead of calling latt_reunitU( struct site *lat ) is that if we already have the gauge transformation matrices malloced we save time by not creating and destroying another temporary
 **/
void
random_transform( struct site *__restrict lat ,
		  GLU_complex *__restrict *__restrict gauge ) ;

/**
   @fn void random_suNC( struct site *__restrict lat )
   @brief generates a random lattice gauge field
   @param lat :: lattice gauge field
 */
void 
random_suNC( struct site *__restrict lat ) ;

/**
   @fn void reunit_gauge_slices( GLU_complex **gauge1 , GLU_complex **gauge2 )
   @brief reunitarise two gauge transformation slices
 */
void
reunit_gauge_slices( GLU_complex **gauge1 , 
		     GLU_complex **gauge2 ) ;

/**
   @fn void trivial( struct site *__restrict lat ) 
   @brief generates a unit gauge field
   @param lat :: lattice gauge field
 */
void 
trivial( struct site *__restrict lat ) ;

#endif
