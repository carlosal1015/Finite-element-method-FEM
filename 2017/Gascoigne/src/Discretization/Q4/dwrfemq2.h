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


#ifndef __DwrFemQ2_h
#define __DwrFemQ2_h

#include "q42d.h"
#include "q43d.h"
#include "baseq2patch.h"
#include "baseq23dpatch.h"
#include "finiteelement.h"
#include "transformation2d.h"
#include "transformation3d.h"

namespace Gascoigne
{

/**********************************************************/

  class DwrFemQ22d : public Q42d
  {
    protected:
      typedef Transformation2d<BaseQ22d>       TransQ2;
      FiniteElement<2,1,TransQ2,BaseQ22dPatch> LowOrderFem;
      HNStructureQ1*                           HNLow;

      void TransformationQ2(FemInterface::Matrix& T, int iq) const;

    public:
      DwrFemQ22d();
      ~DwrFemQ22d();

      void BasicInit(const ParamFile* paramfile);
      void ReInit(const MeshInterface* MP);
  };

/**********************************************************/

  class DwrFemQ23d : public Q43d
  {
    protected:
      typedef Transformation3d<BaseQ23d>       TransQ2;
      FiniteElement<3,2,TransQ2,BaseQ23dPatch> LowOrderFem;
      HNStructureQ1*                           HNLow;

      void TransformationQ2(FemInterface::Matrix& T, int iq) const;

    public:
      DwrFemQ23d();
      ~DwrFemQ23d();

      void BasicInit(const ParamFile* paramfile);
      void ReInit(const MeshInterface* MP);
  };

/**********************************************************/
/**********************************************************/

  class DwrFemQ2Q42d : public DwrFemQ22d
  {
    protected:
      void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex2d& p0, int i, double s) const;

    public:
      DwrFemQ2Q42d() : DwrFemQ22d() { }
      ~DwrFemQ2Q42d() { }

      void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
      void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
      void BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const;

      void MassMatrix(MatrixInterface& M) const;
      void MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const;

      void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

      void HNAverage   (GlobalVector& x)       const { HNLow->Average(x); }
      void HNDistribute(GlobalVector& x)       const { HN->Distribute(x); }
      void HNZero      (GlobalVector& x)       const { HNLow->Zero(x); }
      bool HNZeroCheck (const GlobalVector& x) const { return HNLow->ZeroCheck(x); }
  };

/**********************************************************/

  class DwrFemQ4Q22d : public DwrFemQ22d
  {
    protected:
      void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex2d& p0, int i, double s) const;

    public:
      DwrFemQ4Q22d() : DwrFemQ22d() { }
      ~DwrFemQ4Q22d() { }

      void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
      void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
      void BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const;

      void MassMatrix(MatrixInterface& M) const;
      void MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const;

      void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

      void HNAverage   (GlobalVector& x)       const { HN->Average(x); }
      void HNDistribute(GlobalVector& x)       const { HNLow->Distribute(x); }
      void HNZero      (GlobalVector& x)       const { HN->Zero(x); }
      bool HNZeroCheck (const GlobalVector& x) const { return HN->ZeroCheck(x); }
  };

/**********************************************************/
/**********************************************************/

  class DwrFemQ2Q43d : public DwrFemQ23d
  {
    protected:
      void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex3d& p0, int i, double s) const;

    public:
      DwrFemQ2Q43d() : DwrFemQ23d() { }
      ~DwrFemQ2Q43d() { }

      void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
      void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
      void BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const;

      void MassMatrix(MatrixInterface& M) const;
      void MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const;

      void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

      void HNAverage   (GlobalVector& x)       const { HNLow->Average(x); }
      void HNDistribute(GlobalVector& x)       const { HN->Distribute(x); }
      void HNZero      (GlobalVector& x)       const { HNLow->Zero(x); }
      bool HNZeroCheck (const GlobalVector& x) const { return HNLow->ZeroCheck(x); }
  };

/**********************************************************/

  class DwrFemQ4Q23d : public DwrFemQ23d
  {
    protected:
      void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex3d& p0, int i, double s) const;

    public:
      DwrFemQ4Q23d() : DwrFemQ23d() { }
      ~DwrFemQ4Q23d() { }

      void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
      void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
      void BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const;

      void MassMatrix(MatrixInterface& M) const;
      void MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const;

      void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

      void HNAverage   (GlobalVector& x)       const { HN->Average(x); }
      void HNDistribute(GlobalVector& x)       const { HNLow->Distribute(x); }
      void HNZero      (GlobalVector& x)       const { HN->Zero(x); }
      bool HNZeroCheck (const GlobalVector& x) const { return HN->ZeroCheck(x); }
  };

/**********************************************************/

}

#endif
