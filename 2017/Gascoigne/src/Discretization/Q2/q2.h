/**
*
* Copyright (C) 2004, 2005, 2006, 2009 by the Gascoigne 3D authors
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


#ifndef  __Q2_h
#define  __Q2_h

#include  "patchdiscretization.h"
#include  "hnstructureq1.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q2

////
////
/////////////////////////////////////////////

class Q2 : public PatchDiscretization
{
protected:
  HNStructureQ1*    HN;

  nvector<int> GetLocalIndices(int iq) const {
    return *GetPatchMesh()->IndicesOfPatch(iq);
  }
  void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

public:

//
////  Con(De)structor 
//

  Q2();
  ~Q2();

  int n() const;
  int nc() const;
  int n_withouthanging() const;
  void ReInit(const MeshInterface* MP);

  void Interpolate(GlobalVector& u, const DomainInitialCondition& U) const;
  void StrongDirichletMatrix(MatrixInterface& A, int col, const std::vector<int>& comp) const;
  void StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const std::vector<int>& comp) const;
  void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp, double d) const;
  void StrongDirichletVectorZero(GlobalVector& u, int col, const std::vector<int>& comp) const;
  void InterpolateSolution(GlobalVector& u, const GlobalVector& uold)const;
  void StrongPeriodicVector(GlobalVector& u, const PeriodicData& BF, int col, const std::vector<int>& comp, double d) const;
  void HNAverage   (GlobalVector& x) const;
  void HNDistribute(GlobalVector& x) const;
  void HNZero      (GlobalVector& x) const;
  bool HNZeroCheck (const GlobalVector& x) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const;
  void MassMatrix(MatrixInterface& A) const;
  void Structure(SparseStructureInterface* SI) const;
  void InitFilter(nvector<double>& F) const;
};
}

#endif
