/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (GLU_timer.h) is part of GLU.

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
   @file GLU_timer.h
   @brief function prototypes for the timing of functions
   @warning requires sys/time.h
 */
#ifndef GLU_TIMER_H
#define GLU_TIMER_H

/**
   @fn char* get_date( void )
   @brief returns the current date if time.h exists, otherwise returns nonsense
 **/
char*
get_date( void ) ;

/**
   @fn double print_time( void )
   @brief prints to stdout the time
   Here are the two timing functions I generally use, they are meant to
   be blocked in "#ifdef TIME_GF {...} #endif"

   @return the difference from start_timer()
 **/
double
print_time( void ) ;

/**
   @fn void start_timer( void )
   @brief This refreshes the timer, meant to be followed by a print_time statement.
 **/
void 
start_timer( void ) ;

#endif
