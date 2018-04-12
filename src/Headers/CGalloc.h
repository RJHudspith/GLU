#ifndef GLU_CGALLOC_H
#define GLU_CGALLOC_H

int
allocate_temp_cg( struct gauges *G ,
		  struct CGtemps *CG ,
		  const GLU_bool FACG ) ;

void
free_temp_cg( struct gauges G ,
	      struct CGtemps CG ,
	      const GLU_bool FACG ) ;

int
allocate_temp_lg( struct CGtemps *CG ,
		  const GLU_bool FACG ) ;

void
free_temp_lg( struct CGtemps CG ,
	      const GLU_bool FACG ) ;


#endif
