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


#ifndef  __levelsorter_h
#define  __levelsorter_h

/*---------------------------------------------------*/

namespace Gascoigne
{
class LevelSorter2d
{
protected:

  const HierarchicalMesh2d& HM;

public:

  LevelSorter2d(const HierarchicalMesh2d& HMHM): HM(HMHM) {}
  bool operator()(int i, int j) const
    {
      return ( HM.quad(i).level() > HM.quad(j).level() );
    }
};

/*---------------------------------------------------*/

class HangEdgeSort3
{
protected:

  const LevelMesh2d& LR;

public:

  HangEdgeSort3(const LevelMesh2d& L) : LR(L) {}
  bool operator() (int i, int j) const
    {
      return (!LR.EdgeIsHangingGlobalIndex(i) && LR.EdgeIsHangingGlobalIndex(j));
    }
};
}

#endif
