/**
*
* Copyright (C) 2004, 2006, 2011 by the Gascoigne 3D authors
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


#ifndef  __HNStructureQ22d_h
#define  __HNStructureQ22d_h

#include  "hnstructureq12d.h"

namespace Gascoigne
{

/*-----------------------------------------*/

class HNStructureQ22d : public  HNStructureQ12d
{
protected:
  DoubleVector q1wei;

public:

  HNStructureQ22d();
    
  void Average   (GlobalVector& u) const;
  void Distribute(GlobalVector& u) const;
  void CondenseHanging(EntryMatrix& E, IntVector& indices) const;
  void CondenseHangingLowerHigher(EntryMatrix& E, IntVector& indices) const;
  void CondenseHangingHigherLower(EntryMatrix& E, IntVector& indices) const;
  void CondenseHangingMixed(EntryMatrix& E, IntVector& indices, int k) const;
  void CondenseHanging(IntVector& indices) const;
  void CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const { 
    std::cerr << "\"HNStructureQ22d::CondenseHangingPatch\" not written!" << std::endl;
    abort();
  }

  //void NewCondenseHanging(EntryMatrix& E, IntVector& indices1, IntVector& indices2) const;
};
}
#endif
