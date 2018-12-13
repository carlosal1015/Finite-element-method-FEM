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


#ifndef __hangcontainer2d_h
#define __hangcontainer2d_h

#include  "hanglist.h"
#include  "gascoigne.h"

/*********************************************************************/

namespace Gascoigne
{
class HangContainer2d
{
 protected:

  typedef  fixarray<2,int>               EdgeVector;
  
  HangList<2>   VertexToBeDeleted;  // ehemals haengende Knoten
  HangList<2>   VertexToBeCreated;  // neue haengende Knoten
  HangList<2>   NotAnyMoreHanging;  // Knoten bleibt aber haengt nicht mehr

  HangList<2>&  Hanging;  // dies sind die aktuellen haengenden Knoten

 public:

  HangContainer2d(HangList<2>& lh) : Hanging(lh) {}

  // Zugriff

  int  NToBeDeleted()   const { return VertexToBeDeleted.size(); }
  int  NToBeCreated()   const { return VertexToBeCreated.size(); }

  const HangList<2>& Deleting()   const { return VertexToBeDeleted;}
  const HangList<2>& Creating()   const { return VertexToBeCreated;}
  const HangList<2>& NotAnyMore() const { return NotAnyMoreHanging;}
        HangList<2>& NotAnyMore()       { return NotAnyMoreHanging;}

  bool ToBeDeleted(const EdgeVector& v) const;
  bool ToBeCreated(const EdgeVector& v) const;

  void make_consistent() { VertexToBeCreated.make_consistent(VertexToBeDeleted);}

  void load_elimination(IntVector&) const;

  void update_olds (IntVector&, const IntVector&);
  void update_news (const IntVector&,int);

  int  vertex_index (const EdgeVector&) const;

  void ghost_coarse (EdgeVector&,  int, int);
  void ghost_refine (EdgeVector&,  int);

  void NeighbourSwapper();
};
}

/*********************************************************************/

#endif
