/**
*
* Copyright (C) 2005, 2010, 2011 by the Gascoigne 3D authors
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


/*----------------------------   functionalcontainer.h     ---------------------------*/
/*      $Id$                 */
#ifndef __functionalcontainer_H
#define __functionalcontainer_H
/*----------------------------   functionalcontainer.h     ---------------------------*/


#include <map>
#include <string> 
#include "functional.h"

namespace Gascoigne
{
  
 class  FunctionalContainer : public std::map<std::string, const Functional*>
   {
     public:

     void AddFunctional(const std::string& label, const Functional* P)
       {
	 if (find(label)!=end())
	   {
	     std::cerr << "Functional " << label << " already present!\n";
             abort();
	   }
	 (*this)[label]=P;
       }
     
     void RemoveFunctional(const std::string& label)
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Problemdescriptor " << label << " not present!\n";
             abort();
	   }
	 this->erase(label);
       }
     
     const Functional* GetFunctional(const std::string& label) const
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Functional " << label << " not present!\n";
             abort();
	   }
	 return find(label)->second;
       }

     int GetIndex(const std::string& label) const
     {
       int i = 0;
       for (FunctionalContainer::const_iterator it = this->begin(); it!=this->end();++it,++i)
	 {
	   if(it->first == label)
	     return i;
	 }
       std::cerr<<"Label not found"<<std::endl;
       abort();
     }
     
   };
 
}


/*----------------------------   functionalcontainer.h     ---------------------------*/
/* end of #ifndef __functionalcontainer_H */
#endif
/*----------------------------   functionalcontainer.h     ---------------------------*/
