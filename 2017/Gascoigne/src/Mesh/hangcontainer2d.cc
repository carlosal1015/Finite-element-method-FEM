/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#include "hangcontainer2d.h"
#include "find_in_linehang.h"


using namespace std;

/*********************************************************************/

namespace Gascoigne
{
bool HangContainer2d::ToBeDeleted(const EdgeVector& v) const
{
  return (VertexToBeDeleted.find(v) != VertexToBeDeleted.end());
}

/*********************************************************************/

bool HangContainer2d::ToBeCreated(const EdgeVector& v) const
{
  return (VertexToBeCreated.find(v) != VertexToBeCreated.end());
}

/*********************************************************************/

void HangContainer2d::NeighbourSwapper()
{  
  for(HangList<2>::iterator p=VertexToBeDeleted.begin(); p!=VertexToBeDeleted.end(); p++)
    {
      int r = p->second.rneighbour();
      if (r<0)
	{
	  p->second.rneighbour() = p->second.cneighbour();
	  p->second.cneighbour() = r;
	}
    }
}

/*********************************************************************/

void HangContainer2d::load_elimination(IntVector& v) const
{
  for(HangList<2>::const_iterator p=VertexToBeDeleted.begin(); p!=VertexToBeDeleted.end(); p++)
    {
      v.push_back(p->second.hanging());
    }
}

/*********************************************************************/

void HangContainer2d::update_olds(IntVector& v, const IntVector& c)
{
  VertexToBeCreated.update(v,c);
  VertexToBeDeleted.update(v,c);
  NotAnyMoreHanging.update(v,c);      
  Hanging.update(v,c);      
}

/*********************************************************************/

int HangContainer2d::vertex_index(const EdgeVector& edge) const
{
  HangList<2>::const_iterator p;

  p = VertexToBeCreated.find(edge);
  if(p!=VertexToBeCreated.end()) return p->second.hanging();
  
  p = NotAnyMoreHanging.find(edge);
  if(p!=NotAnyMoreHanging.end()) return p->second.hanging();

  p = Hanging.find(edge);
  if(p!=Hanging.end()) return p->second.hanging();

  return -1;
}

/*********************************************************************/

void HangContainer2d::update_news(const IntVector& vnew, int i)
{
  // cerr << "new_hangs()" << endl;
  // newhangs-hanging setzten fuer new-quad und linehang fuer die zukunft
  
  for (HangList<2>::iterator p = VertexToBeCreated.begin(); p!=VertexToBeCreated.end(); p++)
    {
      p->second.hanging() = vnew[i];
      HangList<2>::iterator Lp = Hanging.find(p->first);
      if(Lp!=Hanging.end())
	{
	  Lp->second.hanging() = vnew[i]; 
	}
      i++;
    }
}

/*********************************************************************/

void HangContainer2d::ghost_coarse(EdgeVector& edge, int f, int edge_vertex)
{
  assert(!find_in_linehang(VertexToBeDeleted,edge).second);

  HangList<2>::iterator e = Hanging.find(edge);
  
  if(e!=Hanging.end())
    {
      /* vertex is hang and has to be deleted */
      if (e->second.cneighbour()==f)
	{
	  swap(e->second.rneighbour(),e->second.cneighbour());
	}
      VertexToBeDeleted.move(Hanging,e);
    }
  else
    {
      Hang  h(edge_vertex,-1,f);
      Hanging.insert(make_pair(edge,h));
    }
}

/*********************************************************************/

void HangContainer2d::ghost_refine(EdgeVector& edge, int f)
{
  HangList<2>::iterator d = VertexToBeDeleted.find(edge);
  HangList<2>::iterator e = Hanging.find(edge);

  if(d!=VertexToBeDeleted.end())
    {
      assert(e==Hanging.end());
      if (d->second.cneighbour()==f)
	{
	  swap(d->second.rneighbour(),d->second.cneighbour());
	}
      Hanging.move(VertexToBeDeleted,d);
    }
  else // not in VertexToBeDeletedang
    {
      if(e!=Hanging.end())
	{
	  NotAnyMoreHanging.move(Hanging,e);
	  // pruefen ob hang in VertexToBeCreatedang
	  HangList<2>::iterator Np = VertexToBeCreated.find(edge);
	  if (Np!=VertexToBeCreated.end())
	    {
	      Np->second.cneighbour() = f;
	    }
	}
      else  // not in VertexToBeDeleted and not Hanging
	{
	  Hang h(-1,f,-1);
	  Hanging.insert(make_pair(edge,h));
	  VertexToBeCreated.insert(make_pair(edge,h));
	}
    }
}
}

/*********************************************************************/
