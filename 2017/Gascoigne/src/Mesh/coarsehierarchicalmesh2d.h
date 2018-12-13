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


#ifndef __coarsehierarchicalmesh2d_h
#define __coarsehierarchicalmesh2d_h

#include  "hierarchicalmesh2d.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class CoarseHierarchicalMesh2d : public HierarchicalMesh2d
{
  protected :
    
  IntSet  CellRefList, CellCoarseList;
  IntVector   cn2o;

  void loop(IntVector& dst);

  public:
  
  CoarseHierarchicalMesh2d(const HierarchicalMesh2d&);
  void BasicInit(int pdepth);
  void GetRefinedList(IntVector&);
  void GetCoarsedList(IntVector&);
  void refine(const IntVector& cell_ref, const IntVector& cell_coarse);
};
}

/*---------------------------------------------------*/

#endif
