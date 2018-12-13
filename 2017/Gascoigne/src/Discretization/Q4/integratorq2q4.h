/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#ifndef __IntegratorQ2Q4_h
#define __IntegratorQ2Q4_h

#include "basicintegrator.h"

namespace Gascoigne
{

/**********************************************************/

  template<int DIM>
  class IntegratorQ2Q4 : public BasicIntegrator
  {
    protected:
      double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}
      int PatchMeshNr2IntegratorNr(int in) const;

    public:
      IntegratorQ2Q4<DIM>() : BasicIntegrator() { }
      ~IntegratorQ2Q4<DIM>() { }

      void BasicInit() { }

      void Form(const Equation& EQ, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U, const LocalData& Q, const LocalData& QC) const;
      void AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U, const LocalData& Q, const LocalData& QC) const;
      void BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U, int ile, int col, LocalData& Q, const LocalData& QC) const;
      void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalData& Q, const LocalData& QC) const;
      void BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, int ile, int col, const LocalData& Q, const LocalData& QC) const;
      void DiracRhsPoint(LocalVector& b, const FemInterface& FemH, const FemInterface& FemL, const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int j, const LocalData& Q, const LocalData& QC) const;

      double MassMatrix(EntryMatrix& E, const FemInterface& FemH, const FemInterface& FemL) const;
      void MassForm(const TimePattern& TP, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U) const;

      std::string GetName() const { return "IntegratorQ2Q4"; }
  };

/**********************************************************/

}

#endif
