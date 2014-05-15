/*
    Copyright 2013 Renwick James Hudspith

    This file (GLU_memcheck.h) is part of GLU.

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
   @file GLU_memcheck.h
   @brief function prototypes to check linux's meminfo to see if we have space for our mallocs

   @warning requires unistd.h so LINUX only for now
 */

#ifndef GLU_MEMCHECK_H
#define GLU_MEMCHECK_H

/**
   @fn short int have_memory_gauge( void )
   @brief memory checking for the gauge field
   checks whether we have enough memory to allocate the gauge field
   checks are turned on with the MEMORY_SAFE_SMEAR define
   @return #GLU_FAILURE or #GLU_SUCCESS
 */
short int
have_memory_gauge( void ) ;

/**
   @fn short int have_memory_hyp( const struct sm_info SMINFO )
   @brief memory checking for our hyp-blocking routines
   @param SMINFO :: general smearing information

   checks to see which Hypercubic blocking routine is suitable
   for the amount of memory we are going to allocate.
   @return fast, medium or slow #GLU_speed
 **/
short int 
have_memory_hyp( const struct sm_info SMINFO ) ;

/**
   @fn short int have_memory_wf( const struct sm_info SMINFO )
   @brief choose which wilson flow routine we should use

   checks to see whether we can use the slightly faster wilson flow
   routine or whether we have to use the slower but memory-cheaper
   version of the code. 
   @return fast or slow #GLU_speed
 **/
short int 
have_memory_wf( const struct sm_info SMINFO ) ;

/**
   @fn short int have_memory_Cgf( )
   @brief Checks to see whether we can perform Coulomb gauge fixing
   @return #GLU_FAILURE or #GLU_SUCCESS
 **/
short int 
have_memory_Cgf( void ) ;

/**
   @fn short int have_memory_Lgf( )
   @brief Checks to see whether we can use the slightly faster Landau gauge
   fixing code.  

   The difference between the two is slight "O(10%)"
   @return fast or slow #GLU_speed
 **/
short int 
have_memory_Lgf( void ) ;

/**
   @fn short int have_memory_readers_writers( const GLU_output config_type ) ;
   @param config_type :: which configuration we are reading or writing.
   @brief Checks to see whether we can use the slightly faster Landau gauge
   fixing code.  

   @return fast or slow #GLU_speed
 **/
short int 
have_memory_readers_writers( const GLU_output config_type ) ;

#endif
