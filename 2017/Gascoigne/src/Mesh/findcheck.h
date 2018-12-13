/**
*
* Copyright (C) 2004, 2011 by the Gascoigne 3D authors
*
* This file is part of Gascoigne 3D
*
* Gascoigne 3D is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* Gascoigne 3D is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT for further information
* on this license.
*
**/


#ifndef  __findcheck_h
#define  __findcheck_h

/*------------------------------------------*/

#define FindMacro(where,message) \
    { CHashIt ip = ##where .find(i); \
       if(ip==##where .end()) { \
          std::string s = "##message ";\
	  cerr << "Index::" << s << " g2l()\n";\
	  cerr << "there is no " << s << " "<< i << std::endl;\
	  abort();}\
      return ip->second; }

#define CheckMacro(where) \
    { CHashIt ip = ##where .find(i);\
      if(ip==##where .end()) return -2;\
      return ip->second;}

/*------------------------------------------*/

#endif
