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


#ifndef __EdgeInfo_h
#define __EdgeInfo_h

#include "edge.h"
#include "gascoigne.h"

/**********************************************************/

namespace Gascoigne
{
template<int DIM>
class EdgeInfo
{
 protected:

  int                    _count;
  fixarray<2*DIM-2,int>  _vertex;
  LocalVector _u;
  const Edge*            _edge;

 public:

  EdgeInfo<DIM>() {}
  ~EdgeInfo<DIM>() {}

  void BasicInit(const Edge*, int, const fixarray<2*DIM-2,int>&);
  void AddNodes(const LocalVector&);

  const fixarray<2*DIM-2,int>&  GetVertex() const { return _vertex; }
  const LocalVector& GetValue()  const { return _u; }
  const Edge&                   GetEdge()   const { return *_edge; }
  int                           GetCount()  const { return _count; }
  fixarray<2*DIM-2,double>      GetNorm()   const;

  void ShowStatistics() const;
};
}

/**********************************************************/

#endif
