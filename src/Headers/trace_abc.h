/**
   @file trace_abc.h
   @brief prototype declarations for the trace of product of 3 matrices
 */
#ifndef TRACE_ABC_H
#define TRACE_ABC_H

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
