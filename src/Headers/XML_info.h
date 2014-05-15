/*
    Copyright 2013 Renwick James Hudspith

    This file (XML_info.h) is part of GLU.

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
   @file XML_info.h
   @brief has the xml information parser for scidac and ILDG xml
 */
#ifndef GLU_XML_INFO_H
#define GLU_XML_INFO_H

/**
   @fn int parse_and_set_xml_SCIDAC( char *xml_info , struct head_data *HEAD_DATA )
   @param xml_info :: the c-string of xml tags for this part of the header
   @param HEAD_DATA :: the struct of header data, some things are written out to

   @warning The lattice geometry in the struct Latt is written in here

   It takes the xml info as a string, and tokenizes the <> tags. Then looks for 
   specific tag names.
 */
int
parse_and_set_xml_SCIDAC( char *xml_info ,
			  struct head_data *HEAD_DATA ) ;

#endif
