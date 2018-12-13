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


#ifndef __levelcomparer2d_h
#define __levelcomparer2d_h

#include "compareclass.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class LevelComparer2d
{
  const HierarchicalMesh2d&   Mesh;
  const IntVector& v;

  public:
  
  LevelComparer2d(const HierarchicalMesh2d& HM, const IntVector& vv) : 
    Mesh(HM), v(vv) {};

  int size() const { return v.size(); }
  int operator[](int i) const 
    {
      return Mesh.quad(v[i]).level();
    }
};
}

/*---------------------------------------------------*/

#endif
