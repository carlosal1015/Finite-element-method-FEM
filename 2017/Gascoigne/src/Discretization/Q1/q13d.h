/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2008, 2011 by the Gascoigne 3D authors
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


#ifndef  __Q1Simple3d_h
#define  __Q1Simple3d_h

#include  "q1.h"
#include  "edgeinfocontainer.h"
#include  "energyestimatorintegrator.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Simple3d
////
////
/////////////////////////////////////////////

class Q13d : public Q1
{
 protected:

  HNStructureInterface* NewHNStructure();
  
  void EEJumps(EdgeInfoContainer<3>& EIC, const GlobalVector& u, const EnergyEstimatorIntegrator<3>& EEI, const HierarchicalMesh3d* HM) const;
  void EEJumpNorm(EdgeInfoContainer<3>& EIC, DoubleVector& eta, const EnergyEstimatorIntegrator<3>& EEI, const HierarchicalMesh3d* HM) const;
  void EEResidual(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const EnergyEstimatorIntegrator<3>& EEI) const;
  int GetCellNumber(const Vertex3d& p0, Vertex3d& p, int c0=0) const;
  void VertexTransformation(const Vertex3d& p0, Vertex3d& p, int iq) const;

public:

  //
  ////  Con(De)structor 
  //
  
  Q13d();

  std::string GetName() const {return "Q13d";}
  
  void BasicInit(const ParamFile* pf);

  void Interpolate(GlobalVector& u, const DomainInitialCondition& U) const;
  void InterpolateSolutionByPatches(GlobalVector& u, const GlobalVector& uold) const;
  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp, double d) const;
  void StrongPeriodicVector(GlobalVector& u, const PeriodicData& BF, int col, const std::vector<int>& comp, double d) const;

  void EnergyEstimator(EdgeInfoContainerInterface& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const std::string & s_energytype, double d_visc) const;

  nmatrix<double> GetLocalInterpolationWeights() const;
};
}

#endif
