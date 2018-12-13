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


#ifndef  __HNStructureQ13d_h
#define  __HNStructureQ13d_h

#include  "hnstructureq12d.h"
#include  "entrymatrix.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class HNStructureQ13d : public HNStructureQ12d
{
protected:

  typedef   fixarray<9,int>    FaceVector;

  typedef   std::map<int,FaceVector>::iterator        fiterator;
  typedef   std::map<int,FaceVector>::const_iterator  const_fiterator;

  const std::map<int,FaceVector>*        faces;

  fixarray<12,fixarray<3,int> >  lnoe;
  fixarray< 6,fixarray<5,int> >  lnop;

  void CondenseHanging2er(IntVector& indices) const;
  void CondenseHanging4er(IntVector& indices) const;

  void CondenseHanging2er(EntryMatrix& E, IntVector& indices) const;
  void CondenseHanging4er(EntryMatrix& E, IntVector& indices) const;

  fixarray<4,int> GetHangingFace(int i) const;
  fixarray<2,int> GetHangingEdge(int i) const;

public:

  ~HNStructureQ13d();
  HNStructureQ13d();

  int nhnodes() const { return edges->size() + faces->size();} 

  void ReInit(const MeshInterface* m);
  int   hanging(int i) const;

  void MatrixDiag(int ncomp, MatrixInterface& A) const;
  void SparseStructureDiag(SparseStructure* A) const;
  
  void Average(GlobalVector& u) const;
  void Distribute(GlobalVector& u) const;
  void Zero(GlobalVector& u) const;
  bool ZeroCheck(const GlobalVector& u) const;
  
  void Couplings(IntVector& indices) const;
  
  void CondenseHanging(IntVector& indices) const;
  void CondenseHanging(EntryMatrix&, IntVector&) const;
  void CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const;
};
}

#endif
