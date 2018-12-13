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


#ifndef  __HNStructureQ12d_h
#define  __HNStructureQ12d_h

#include  "hnstructureq1.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class HNStructureQ12d : public HNStructureQ1
{
protected:

  const std::map<int,EdgeVector>*        edges;
  DoubleVector                        wei;
  fixarray<4,fixarray<3,int> >           lnoe, lnop;

  double weight(int i) const { return wei[i];}
  int hanging(int i) const;
  const EdgeVector& regular_nodes(int i) const;
  
public:

  HNStructureQ12d();
  ~HNStructureQ12d() {}
  int nhnodes() const {return edges->size();} 
  void SparseStructureDiag(SparseStructure* S) const;
  void ReInit(const MeshInterface* m);
  
  void MatrixDiag(int ncomp, MatrixInterface& A) const;
  void Average(GlobalVector& u) const;
  void Distribute(GlobalVector& u) const;
  void Zero(GlobalVector& u) const;
  bool ZeroCheck(const GlobalVector& u) const;
  
  void CondenseHanging(IntVector& indices) const;
  void CondenseHanging(EntryMatrix& E, IntVector& indices) const;
  void CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const;
};
}

#endif
