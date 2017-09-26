/*
    Copyright 2013-2017 Renwick James Hudspith

    This file (Qmoments.h) is part of GLU.

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
   @file Qmoments.h
   @biref prototype declarations for Qmoment routines
 */
#ifndef GLU_QMOMENTS_H
#define GLU_QMOMENTS_H

/**
   @fn int compute_Qmoments( struct site *__restrict lat , struct Qmoments *Qmom , const struct cut_info CUTINFO , const size_t measurement )
   @brief compute moments of the topological charge
 */
int
compute_Qmoments( struct site *__restrict lat ,
		  struct Qmoments *Qmom , 
		  const struct cut_info CUTINFO ,
		  const size_t measurement ) ;

#endif
