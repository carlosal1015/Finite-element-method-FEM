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


#include  "index.h"


using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne
{
Index::Index() {}

/*---------------------------------------------------*/

Index::Index(const Index& I) 
{
  *this = I;
}

/*---------------------------------------------------*/

Index& Index::operator=(const Index& I)
{
  vl2g  = I.Vertexl2g();
  vg2l  = I.Vertexg2l() ;
  el2g  = I.Edgel2g();
  eg2l  = I.Edgeg2l();
  hl2g  = I.Hexl2g();
  hg2l  = I.Hexg2l();
  ql2g  = I.Quadl2g();
  qg2l  = I.Quadg2l();

  return *this;
}

/*---------------------------------------------------*/

ostream& operator<<(ostream& os, const Index& I)
{
  os << "Vertex l2g " << I.VertexSize()<<endl;
  os << I.Vertexl2g();
  os << "Vertex g2l " << I.VertexSize()<<endl;
  os << I.Vertexg2l() << " ";

  os << "Edge l2g " << I.EdgeSize()<<endl;
  os << I.Edgel2g();
  os << "Edge g2l " << I.EdgeSize()<<endl;
  os << I.Edgeg2l() << " ";

  os << "Hex l2g " << I.HexSize()<<endl;
  os << I.Hexl2g();
  os << "Hex g2l " << I.HexSize()<<endl;
  os << I.Hexg2l() << " ";

  os << "Quad l2g " << I.QuadSize()<<endl;
  os << I.Quadl2g();
  os << "Quad g2l " << I.QuadSize()<<endl;
  os << I.Quadg2l() << " ";

  return os;
}

/*---------------------------------------------------*/

void Index::InitNodes(const IntSet& nodes)
{
  typedef IntSet::const_iterator  iterator;
  int n = nodes.size();
  vl2g.memory(n);
  int count=0;
  iterator p=nodes.begin();
  for(int i=0;i<n;i++)
    {
      vl2g[count] = *p++;
      count++;
    }
  vg2l.clear();
  for(int i=0;i<n;i++)
    {
      vg2l.insert(make_pair(vl2g[i],i));
    }
}

/*---------------------------------------------------*/

void Index::InitEdges(const IntSet& edges)
{
  typedef IntSet::const_iterator  iterator;
  int  n = edges.size();
  el2g.memory(n);
  int count=0;
  for(iterator p=edges.begin();p!=edges.end();p++)
    {
      int i = *p;
      el2g[count++] = i;
    }
}

/*---------------------------------------------------*/

void Index::InitQuads()
{
  int  n  = ql2g.size();
  qg2l.clear();
  for(int i=0;i<n;i++)
    {
      qg2l.insert(make_pair(ql2g[i],i));
    }
}

/*---------------------------------------------------*/

void Index::InitHexs()
{
  int  n  = hl2g.size();
  hg2l.clear();
  for(int i=0;i<n;i++)
    {
      hg2l.insert(make_pair(hl2g[i],i));
    }
}
}

/*---------------------------------------------------*/
