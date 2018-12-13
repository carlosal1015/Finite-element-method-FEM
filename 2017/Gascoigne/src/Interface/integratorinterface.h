/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2010, 2011 by the Gascoigne 3D authors
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


#ifndef  __IntegratorInterface_h
#define  __IntegratorInterface_h


#include  "gascoigne.h"
#include  "equation.h"
#include  "exactsolution.h"
#include  "feminterface.h"
#include  "entrymatrix.h"
#include  "domainrighthandside.h"
#include  "domainfunctional.h"
#include  "boundaryfunctional.h"
#include  "boundaryrighthandside.h"
#include  "boundaryequation.h"
#include  "diracrighthandside.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments IntegratorInterface

  ///
  ///
  /////////////////////////////////////////////

  class IntegratorInterface
  {
    private:

    protected:

    public:
      IntegratorInterface() {}
      virtual ~IntegratorInterface() {}
  
      virtual std::string GetName() const=0;

      virtual void BasicInit() {
        std::cerr << "\"IntegratorInterface::BasicInit\" not written!" << std::endl;
        abort();
      }
      
      virtual void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, 
          const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::Rhs\" not written!" << std::endl;
        abort();
      }
      virtual void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
          const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::Form\" not written!" << std::endl;
        abort();
      }
      virtual void AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FEM, 
          const LocalVector& U, const LocalData& Q, const LocalData& QC) const { 
        std::cerr << "\"IntegratorInterface::AdjointForm\" not written!" << std::endl;
        abort();
      }

      virtual void BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, int ile, 
          int col, const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::BoundaryRhs\" not written!" << std::endl;
        abort();
      }
      virtual void BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
          int ile, int col, const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::BoundaryForm\" not written!" << std::endl;
        abort();
      }
      virtual void BoundaryMatrix (const BoundaryEquation& BE, EntryMatrix& E, const FemInterface& FEM, 
          const LocalVector& U, int ile, int col, const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::BoundaryMatrix\" not written!" << std::endl;
        abort();
      }
      virtual void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
          const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::Matrix\" not written!" << std::endl;
        abort();
      }
      virtual double MassMatrix(EntryMatrix& E, const FemInterface& FEM) const {
        std::cerr << "\"IntegratorInterface::MassMatrix\" not written!" << std::endl;
        abort();
      } 
      virtual void BoundaryMassMatrix (EntryMatrix& E, const FemInterface& FEM, int ile) const
      {
        std::cerr << "\"IntegratorInterface::BoundaryMassMatrix\" not written!" << std::endl;
        abort();
      }

      virtual void MassForm(const TimePattern& TP, LocalVector& F, const FemInterface& FEM, const LocalVector& U) const {
        std::cerr << "\"IntegratorInterface::MassForm\" not written!" << std::endl;
        abort();
      }

      virtual double ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, 
          const LocalVector& U, const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::ComputeDomainFunctional\" not written!" << std::endl;
        abort();
      }
      virtual double ComputeErrorDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, 
          const LocalVector& U, const LocalData& Q, const LocalData& QC) const {
        return ComputeDomainFunctional(F, FEM, U, Q, QC);
      }
      virtual double ComputeBoundaryFunctional(const BoundaryFunctional& F, const FemInterface& FEM, int ile,
          int col, const LocalVector& U, const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::ComputeBoundaryFunctional\" not written!" << std::endl;
        abort();
      }

        
      virtual void EvaluateCellRightHandSide(LocalVector &F, const DomainRightHandSide& CF,const FemInterface& FEM, 
          const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::EvaluateCellRightHandSide\" not written!" << std::endl;
        abort(); }

      virtual void EvaluateBoundaryCellRightHandSide(LocalVector& F, const BoundaryRightHandSide& CF,const FemInterface& FEM, int ile, int col, 
	   const LocalData& Q, const LocalData& QC) const{
        std::cerr << "\"IntegratorInterface::EvaluateBoundaryCellRightHandSide\" not written!" << std::endl;
        abort(); }

      virtual void DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex2d& p, const DiracRightHandSide& DRHS, 
          int i, const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::DiracRhsPoint\" not written!" << std::endl;
        abort();
      }
      virtual void DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex3d& p, const DiracRightHandSide& DRHS, 
          int i, const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::DiracRhsPoint\" not written!" << std::endl;
        abort();
      }
      virtual double ComputePointValue(const FemInterface& E, const Vertex2d& p, const LocalVector& U, int comp) const {
        std::cerr << "\"IntegratorInterface::ComputePointValue\" not written!" << std::endl;
        abort();
      }
      virtual double ComputePointValue(const FemInterface& E, const Vertex3d& p, const LocalVector& U, int comp) const {
        std::cerr << "\"IntegratorInterface::ComputePointValue\" not written!" << std::endl;
        abort();
      }
      virtual void ErrorsByExactSolution(LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, 
          const LocalVector& U, const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::ErrorsByExactSolution\" not written!" << std::endl;
        abort();
      }
      virtual void IntegrateMassDiag(DoubleVector& F, const FemInterface& FEM) const {
        std::cerr << "\"IntegratorInterface::IntegrateMassDiag\" not written!" << std::endl;
        abort();
      }

      virtual void IntegrateBoundaryMassDiag(DoubleVector& F, const FemInterface& FEM, int ile, int col) const {
        std::cerr << "\"IntegratorInterface::IntegrateBoundaryMassDiag\" not written!" << std::endl;
        abort();
      }

      virtual void RhsCurve(LocalVector& F, const FemInterface& FEM, Vertex2d& xr0, Vertex2d& xr1, double H, double ND0,double ND1, int ncomp, int comp)const {
        std::cerr << "\"IntegratorInterface::RhsCurve\" not written!" << std::endl;
        abort();
      }

      virtual void MassMatrixVector(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
          const LocalData& Q, const LocalData& QC) const {
        std::cerr << "\"IntegratorInterface::MassMatrixVector\" not written!" << std::endl;
						assert(0);
      }

  };
}

#endif
