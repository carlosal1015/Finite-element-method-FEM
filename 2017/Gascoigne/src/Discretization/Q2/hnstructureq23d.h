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


#ifndef  __HNStructureQ23d_h
#define  __HNStructureQ23d_h

#include  "hnstructureq13d.h"

namespace Gascoigne
{

/*-----------------------------------------*/

class HNStructureQ23d : public  HNStructureQ13d
{
  fixarray<9,double>  fwei, fq1wei;
  DoubleVector        q1wei;

  fixarray<12,fixarray<3,int> >  lnoe;
  fixarray< 6,fixarray<5,int> >  lnop;

  void CondenseHanging2er(EntryMatrix& E, nvector<int>& indices) const;
  void CondenseHanging2erLowerHigher(EntryMatrix& E, nvector<int>& indices) const;
  void CondenseHanging2erHigherLower(EntryMatrix& E, nvector<int>& indices) const;
  void CondenseHanging4er(EntryMatrix& E, nvector<int>& indices) const;
  void CondenseHanging4erLowerHigher(EntryMatrix& E, nvector<int>& indices) const;
  void CondenseHanging4erHigherLower(EntryMatrix& E, nvector<int>& indices) const;

public:

  HNStructureQ23d();

  void Average(GlobalVector& u) const;
  void Distribute(GlobalVector& u) const;
  void CondenseHanging(EntryMatrix& E, IntVector& indices) const;
  void CondenseHangingLowerHigher(EntryMatrix& E, IntVector& indices) const;
  void CondenseHangingHigherLower(EntryMatrix& E, IntVector& indices) const;
  void CondenseHanging(IntVector& indices) const;
  void CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const {
    std::cerr << "\"HNStructureQ23d::CondenseHangingPatch\" not written!" << std::endl;
    abort();
  }
};

}
#endif
