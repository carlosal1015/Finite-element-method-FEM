/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#include "ghostvectoragent.h"
#include "stlio.h"

namespace Gascoigne
{

/*-------------------------------------------------*/

GhostVectorAgent::GhostVectorAgent() {}

/*-------------------------------------------------*/

GhostVectorAgent::~GhostVectorAgent()
{
  for (iterator p=begin(); p!=end(); p++)
    { 
      if(p->second) 
	{
	  //  Loesche p->first
	  delete p->second; 
	  p->second = NULL;
	} 
    }
}
  
/*-------------------------------------------------*/

void GhostVectorAgent::Register(const VectorInterface& mg) 
{
  iterator p = find(mg);
  if(p==end())
    {
      insert(std::make_pair(mg,static_cast<GlobalVector*>(NULL)));
    }
}
  
/*-------------------------------------------------*/

void GhostVectorAgent::Delete(VectorInterface& mg) 
{
  iterator p=find(mg);
  if (p!=end())
    {
      delete p->second; 
      erase(p);
    }
}
    
/*-------------------------------------------------*/
  
GlobalVector& GhostVectorAgent::operator()(const VectorInterface& g) 
{
  iterator p = find(g);
  if (p==end())
    {
      std::cerr << __FILE__ << ":" << __LINE__;
      std::cerr << ": GhostVectorAgent::operator(): ERROR"<<std::endl;
      std::cerr << __FILE__ << ":" << __LINE__;
      std::cerr << ": Ghostvector '"<< g <<"' not found in list of: "<<std::endl;
      std::cerr << " "<< *this << std::endl;
      abort();
    }
  GlobalVector* vp = p->second;
  if (vp==NULL) 
    {
      std::cerr <<  "GhostVectorAgent  GlobalVector* NULL\t" << p->first;
      std::cerr << "\n" << *this << std::endl;
      abort();
    }
  return *vp;
}

}
