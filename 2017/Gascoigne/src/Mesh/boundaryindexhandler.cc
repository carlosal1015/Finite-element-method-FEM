/**
*
* Copyright (C) 2004, 2005, 2007, 2009, 2011 by the Gascoigne 3D authors
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


#include  "boundaryindexhandler.h"
#include  "set2vec.h"
#include  <set>
#include  "stlio.h"

using namespace std;

/*--------------------------------------------------------------*/

namespace Gascoigne
{
ostream& operator<<(ostream &s, const BoundaryIndexHandler& A)
{
  const IntSet& colors = A.GetColors();
  cerr << "All Colors: " <<  colors << endl;
  for(IntSet::const_iterator p=colors.begin();p!=colors.end();++p)
    {
      cerr << "color: " << *p;
      cerr << "\n\tVertices: " << A.Verteces(*p);
      cerr << "\n\tCells: " << A.Cells(*p);
      cerr << "\n\tLocalind: " << A.Localind(*p);
      cerr << endl;
    }
  return s;
}

/*--------------------------------------------------------------*/

void BoundaryIndexHandler::check() const
{
  const IntSet& colors = GetColors();
  cerr << "All Colors: " <<  colors << endl;
  for(IntSet::const_iterator p=colors.begin();p!=colors.end();++p)
    {
      cerr << "color: " << *p;
      const IntVector& v= Verteces(*p);
      for(int i=0;i<v.size();i++)
	{
	  if(v[i]<0)
	    {
	      cerr << "....BoundaryIndexHandler::check() ERROR\n";
              abort();
	    }
	}
    }
  cerr << endl;
}

/*--------------------------------------------------------------*/

void BoundaryIndexHandler::Equal
(const IntSet& col, const VecMap& v, const VecMap& c, const VecMap& l)
{
  AllColors = col;
  verteces  = v;
  cells     = c;
  localci   = l;
}


/*--------------------------------------------------------------*/

void BoundaryIndexHandler::CopySetToVector
(const vector<IntSet>& H, const IntVector& colorvec, VecMap& dst) const
{
  for (int i=0; i<H.size(); i++)
    {
      IntVector v;
      Set2Vec(v,H[i]);
      int color = colorvec[i];
      dst.insert(make_pair(color,v));
    }
}

/*--------------------------------------------------------------*/

void BoundaryIndexHandler::clear()
{
  AllColors.clear();
  verteces.clear();
  cells   .clear();
  localci .clear();
}

/*--------------------------------------------------------------*/


const IntVector& BoundaryIndexHandler::Verteces(int color) const 
{
  VecMap::const_iterator p = verteces.find(color);
  if(p==verteces.end())
    {
      cerr << "BoundaryIndexHandler::Vertices\tcolor not found: "
	   << color << endl;
      abort();
    }
  return p->second;
}

/*--------------------------------------------------------------*/

const IntVector& BoundaryIndexHandler::Cells(int color) const 
{ 
  VecMap::const_iterator p = cells.find(color);
  if(p==cells.end())
    {
      cerr << "BoundaryIndexHandler::Cells\tcolor not found: "<<color << endl;
      abort();
    }
  return p->second;
}

/*--------------------------------------------------------------*/

const IntVector& BoundaryIndexHandler::Localind(int color) const 
{ 
  VecMap::const_iterator p = localci.find(color);
  if(p==localci.end())
    {
      cerr << "BoundaryIndexHandler::Localind\tcolor not found: "<<color << endl;
      abort();
    }
  return p->second;
}

/*--------------------------------------------------------------*/

const IntVector& BoundaryIndexHandler::Patches(int color) const
{
  VecMap::const_iterator p = patches.find(color);
  if(p==patches.end())
    {
      cerr << "BoundaryIndexHandler::Patches\tcolor not found: "<<color << endl;
      abort();
    }
  return p->second;
}

/*--------------------------------------------------------------*/
  
const IntVector& BoundaryIndexHandler::LocalPatchind(int color) const
{
  VecMap::const_iterator p = localpi.find(color);
  if(p==localpi.end())
    {
      cerr << "BoundaryIndexHandler::LocalPatchind\tcolor not found: "<<color << endl;
      abort();
    }
  return p->second;
}

/*--------------------------------------------------------------*/

void BoundaryIndexHandler::SetPeriodicPairs(std::map<int,std::map<int,int> > mm_PeriodicPairs)
{
  _PeriodicPairs = mm_PeriodicPairs;
}

/*--------------------------------------------------------------*/

const std::map<int,std::map<int,int> > BoundaryIndexHandler::GetPeriodicPairs() const
{
  return _PeriodicPairs;
}

}
