/**
   @file KPHB.h
   @brief prototype functions for quenched heat bath updates
 */
#ifndef KPHB_H
#define KPHB_H

/**
 */
int
hb_update( struct site *lat ,
	   const struct hb_info HBINFO ,
	   const char *traj_name ,
	   const GLU_output storage , 
	   const char *output_details ,
	   const GLU_bool continuation ) ;

#endif
