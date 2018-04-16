/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (GLU_timer.c) is part of GLU.

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
   @file GLU_timer.c
   @brief timing functions
 */

#include "Mainfile.h"

#ifdef HAVE_TIME_H
#include <time.h>
#endif

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

// functions for the elapsed time of a process
#ifdef HAVE_SYS_TIME_H

/**
   @struct GLUtimer
   @brief the timer from sys/time.h
   
   @var t1
   @brief the time in microseconds in double precision
 */
static struct timeval GLUtimer ;
static double t1 ;

// print to stdout the time elapsed in seconds, minutes or hours
double
print_time( void )
{
  gettimeofday( &GLUtimer , NULL ) ;
  double diff = GLUtimer.tv_sec + ( GLUtimer.tv_usec/ 1E6 ) - t1 ;
  fprintf( stdout , "\n[TIMER] elapsed :: %e s \n", diff ) ;
  return diff ;
}

// start the timer
void
start_timer( void )
{
  gettimeofday( &GLUtimer , NULL ) ;
  t1 = GLUtimer.tv_sec + ( GLUtimer.tv_usec / 1E6 ) ;
  return ;
}
#else // do nothing

double
print_time( void ) { return 0.0 ; }

void
start_timer( void ) { return ; }

#endif

// functions for the date
#ifdef HAVE_TIME_H

// gives us GMT, why? Because I am british
char*
get_date( void )
{
  // Puts a newline at the end of the char! Why?
  time_t today = time( NULL ) ;
  char *test2 = asctime( gmtime( &today ) ) ;
  test2[ strlen( test2 ) - 1 ] = '\0' ; // dirty strip the trailing newline
  return test2 ;
}

#else

// seems like a reasonable date
char*
get_date( void )
{
  return "Mon Jan  0 00:00:00 1970" ;
}

#endif
