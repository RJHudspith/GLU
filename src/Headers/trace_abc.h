/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (trace_abc.h) is part of GLU.

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
   @file trace_abc.h
   @brief prototype declarations for the trace of product of 3 matrices
 */
#ifndef GLU_TRACE_ABC_H
#define GLU_TRACE_ABC_H

/**
   @fn double Re_trace_abc_dag_suNC( const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] ) 
   @brief real part of the trace of the product of 3 SU(N) matrices
 */
double
Re_trace_abc_dag_suNC( const GLU_complex a[ NCNC ] , 
		       const GLU_complex b[ NCNC ] , 
		       const GLU_complex c[ NCNC ] ) ;

/**
   @fn void trace_abc( GLU_complex *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief trace of the product of 3 #NCNC matrices
 */
void
trace_abc( GLU_complex *__restrict tr , 
	   const GLU_complex a[ NCNC ] , 
	   const GLU_complex b[ NCNC ] , 
	   const GLU_complex c[ NCNC ] ) ;

/**
   @fn void trace_abc_dag( GLU_complex *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief trace of abc where c is daggered
 */
void
trace_abc_dag( GLU_complex *__restrict tr , 
	       const GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC ] , 
	       const GLU_complex c[ NCNC ] ) ;

/**
   @fn void trace_abc_dag_Re( GLU_real *__restrict tr , const GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief real part of the trace of abc where c is daggered
 */
void
trace_abc_dag_Re( GLU_real *__restrict tr , 
		  const GLU_complex a[ NCNC ] , 
		  const GLU_complex b[ NCNC ] , 
		  const GLU_complex c[ NCNC ] ) ;

#endif
