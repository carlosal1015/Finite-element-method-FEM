/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#ifndef  __Q1_h
#define  __Q1_h


#include  "celldiscretization.h"
#include  "hnstructureinterface.h"
#include  "edgeinfocontainer.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1

////
////
/////////////////////////////////////////////

class Q1 : public CellDiscretization
{
protected:

  HNStructureInterface*    HN;

  IntVector GetLocalIndices(int iq) const {
    return GetMesh()->IndicesOfCell(iq);
  }
  void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;
  
  virtual HNStructureInterface* NewHNStructure()=0;

public:

//
////  Con(De)structor 
//

  Q1();
  ~Q1();

  int n() const {
    return GetMesh()->nnodes();
  }
  int nc() const {
    return GetMesh()->ncells();
  }
  //  const HNStructureQ1* GetHNStructure() const { return HN;}

  void ReInit   (const MeshInterface* MP);
  void StrongDirichletMatrix       (MatrixInterface& A, int col, const std::vector<int>& comp) const;
  void StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const std::vector<int>& comp) const;
  void StrongDirichletVectorZero(GlobalVector& u, int col, const std::vector<int>& comp) const;
  void InterpolateSolution(GlobalVector& u, const GlobalVector& uold)const;
  virtual void InterpolateSolutionByPatches(GlobalVector& u, const GlobalVector& uold) const {}
  void HNAverage   (GlobalVector& x) const;
  void HNDistribute(GlobalVector& x) const;
  void HNZero      (GlobalVector& x) const;
  bool HNZeroCheck (const GlobalVector& x) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const;
  void MassMatrix(MatrixInterface& A) const;
  void Structure(SparseStructureInterface* SI) const;
  void InitFilter(DoubleVector& F) const;
  virtual void EnergyEstimator(EdgeInfoContainerInterface& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const std::string & s_energytype, double d_visc) const=0;
};
}


#endif
