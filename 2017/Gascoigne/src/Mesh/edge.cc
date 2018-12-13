/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#include "edge.h"

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne
{
Edge& Edge::operator=(const Edge& e)
{
  c1 = e.master();
  c2 = e.slave();
  l1 = e.LocalMasterIndex();
  l2 = e.LocalSlaveIndex ();

  return *this;
}

/*---------------------------------------------------*/

pair<int,int> Edge::EdgeNeighbour(int i) const
{
  //const Edge&  E = edge(quad(i).edge(e));
  int in = master();
  int il = LocalMasterIndex();
  if(in==i)    
    {
      in = slave();
      il = LocalSlaveIndex();
    }
  return make_pair(in,il);
}

/*---------------------------------------------------*/

void Edge::swapping(int newindex)
{
  master() = newindex;
  LocalMasterIndex() = LocalSlaveIndex();
  LocalSlaveIndex()  = -1;
  slave() = -1;
}

/*---------------------------------------------------*/

void Edge::setmaster(int newindex, int newlocal)
{
  master() = newindex;
  LocalMasterIndex() = newlocal;
  LocalSlaveIndex()  = -1;
  slave() = -1;
}

/*---------------------------------------------------*/

void Edge::BinWrite(ostream &s) const
{
  int sizeInt = sizeof(int);
  s.write(reinterpret_cast<const char*>(&c1),sizeInt);
  s.write(reinterpret_cast<const char*>(&l1),sizeInt);
  s.write(reinterpret_cast<const char*>(&c2),sizeInt);
  s.write(reinterpret_cast<const char*>(&l2),sizeInt);
}

/*---------------------------------------------------*/

void Edge::BinRead(istream &s)
{
  int sizeInt = sizeof(int);
  s.read(reinterpret_cast<char*>(&c1),sizeInt);
  s.read(reinterpret_cast<char*>(&l1),sizeInt);
  s.read(reinterpret_cast<char*>(&c2),sizeInt);
  s.read(reinterpret_cast<char*>(&l2),sizeInt);
}

/*---------------------------------------------------*/

ostream& operator<<(ostream &s, const Edge& A)
{
  s << A.master()  << " ";
  s << A.LocalMasterIndex()  << " ";
  s << A.slave()   << " ";
  s << A.LocalSlaveIndex() << " "<< endl;
  
  return s;
}

/*---------------------------------------------------*/

istream& operator>>(istream &s, Edge& A)
{
  s >> A.master() >> A.LocalMasterIndex() >> A.slave() >> A.LocalSlaveIndex();

  return s;
}
}

/*---------------------------------------------------*/
