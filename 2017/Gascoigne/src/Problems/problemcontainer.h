/**
*
* Copyright (C) 2005, 2011 by the Gascoigne 3D authors
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


/*----------------------------   problemcontainer.h     ---------------------------*/
/*      $Id$                 */
#ifndef __problemcontainer_H
#define __problemcontainer_H
/*----------------------------   problemcontainer.h     ---------------------------*/


#include <map>
#include <string> 
#include "problemdescriptorinterface.h"

namespace Gascoigne
{
  
 class  ProblemContainer : public std::map<std::string, const ProblemDescriptorInterface*>
   {
     public:

     void AddProblem(const std::string& label, const ProblemDescriptorInterface* P)
       {
	 if (find(label)!=end())
	   {
	     std::cerr << "Problemdescriptor " << label << " already present!\n";
             abort();
	   }
	 (*this)[label]=P;
       }
     
     void RemoveProblem(const std::string& label)
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Problemdescriptor " << label << " not present!\n";
             abort();
	   }
	 this->erase(label);
       }
     
     const ProblemDescriptorInterface* GetProblem(const std::string& label) const
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Problemdescriptor " << label << " not present!\n";
             abort();
	   }
	 return find(label)->second;
       }
     
   };
 
}



/*----------------------------   problemcontainer.h     ---------------------------*/
/* end of #ifndef __problemcontainer_H */
#endif
/*----------------------------   problemcontainer.h     ---------------------------*/
