/**
*
* Copyright (C) 2004, 2005, 2006, 2007 by the Gascoigne 3D authors
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


#ifndef  __GalerkinIntegrator_h
#define  __GalerkinIntegrator_h

#include  "basicintegrator.h"
#include  "integrationformula.h"
#include  "integrationformulasummed.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GalerkinIntegrator
////
////
/////////////////////////////////////////////

template<int DIM>
class GalerkinIntegrator : public BasicIntegrator
{
private:

  IntegrationFormulaInterface*  IFF;
  IntegrationFormulaInterface*  IFE;
  IntegrationFormulaInterface*  IFB;
  IntegrationFormulaInterface*  IFM;

protected:

  IntegrationFormulaInterface*& FormFormulaPointer() { return IFF;}
  IntegrationFormulaInterface*& ErrorFormulaPointer() { return IFE;}
  IntegrationFormulaInterface*& MassFormulaPointer() { return IFM;}
  IntegrationFormulaInterface*& BoundaryFormulaPointer() { return IFB;}

  const IntegrationFormulaInterface* FormFormula() const { assert(GalerkinIntegrator<DIM>::IFF); return GalerkinIntegrator<DIM>::IFF;}
  const IntegrationFormulaInterface* MassFormula() const { assert(GalerkinIntegrator<DIM>::IFM); return GalerkinIntegrator<DIM>::IFM;}
  const IntegrationFormulaInterface* ErrorFormula() const { assert(GalerkinIntegrator<DIM>::IFE); return GalerkinIntegrator<DIM>::IFE;}
  const IntegrationFormulaInterface* BoundaryFormula() const { assert(GalerkinIntegrator<DIM>::IFB); return GalerkinIntegrator<DIM>::IFB;}

  double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}

public:

//
////  Con(De)structor 
//

  GalerkinIntegrator<DIM>();
  ~GalerkinIntegrator<DIM>();
  
  std::string GetName() const {return "Galerkin";}
  void BasicInit();

  void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, 
      const LocalData& Q, const LocalData& QC) const;
  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
      const LocalData& Q, const LocalData& QC) const;
  void AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
      const LocalData& Q, const LocalData& QC) const;
  void BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FEM, const LocalVector& U, int ile, int col, 
      const LocalData& Q, const LocalData& QC) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
      const LocalData& Q, const LocalData& QC) const;
  void BoundaryMatrix (const BoundaryEquation& BE, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, int ile, int col, 
      const LocalData& Q, const LocalData& QC) const;
  double MassMatrix(EntryMatrix& E, const FemInterface& FEM) const;
  void BoundaryMassMatrix (EntryMatrix& E, const FemInterface& FEM, int ile) const;
  void MassForm(const TimePattern& TP, LocalVector& F, const FemInterface& FEM, const LocalVector& U) const;

  void RhsPoint(LocalVector& b, const FemInterface& E, const Vertex<DIM>& p, int comp) const;
  void DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int j,
      const LocalData& Q, const LocalData& QC) const;
  double ComputePointValue(const FemInterface& E, const Vertex<DIM>& p, const LocalVector& U, int comp) const;
  double ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U,
      const LocalData& Q, const LocalData& QC) const;
  double ComputeErrorDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U,
      const LocalData& Q, const LocalData& QC) const;
  double ComputeBoundaryFunctional(const BoundaryFunctional& F, const FemInterface& FEM, int ile, int col, const LocalVector& U,
      const LocalData& Q, const LocalData& QC) const;
  void EvaluateCellRightHandSide(LocalVector& F, const DomainRightHandSide& CF,const FemInterface& FEM, 
      const LocalData& Q, const LocalData& QC) const;
  void EvaluateBoundaryCellRightHandSide(LocalVector& F, const BoundaryRightHandSide& CF,const FemInterface& FEM, int ile, int col, 
      const LocalData& Q, const LocalData& QC) const;

  void ErrorsByExactSolution(LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, const LocalVector& U, 
      const LocalData& Q, const LocalData& QC) const;

  void BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, int ile, int col,
      const LocalData& Q, const LocalData& QC) const;

  void IntegrateMassDiag(DoubleVector& F, const FemInterface& FEM) const ;

  void IntegrateBoundaryMassDiag(DoubleVector& F, const FemInterface& FEM, int ile, int col) const ;

  void RhsCurve(LocalVector& F, const FemInterface& FEM, Vertex<DIM>& xr0,
      Vertex<DIM>& xr1, double H, double ND0,double ND1, int ncomp, int comp) const;


};
}


#endif
