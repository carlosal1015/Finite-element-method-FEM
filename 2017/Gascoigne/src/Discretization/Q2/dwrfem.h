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


#ifndef __Dwrfem_h
#define __Dwrfem_h

#include "q22d.h"
#include "q23d.h"
#include "baseq1patch.h"
#include "baseq13dpatch.h"
#include "transformation2d.h"
#include "transformation3d.h"
#include "finiteelement.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class DwrFem2d : public Q22d
{
 protected:
  typedef Transformation2d<BaseQ12d>         TransQ1;
  FiniteElement<2,1,TransQ1,BaseQ12dPatch>   LowOrderFem;
  HNStructureQ1*                             HNLow;
  
  void TransformationQ1(FemInterface::Matrix& T, int iq) const;

 public:
  DwrFem2d();
  ~DwrFem2d();

  void BasicInit(const ParamFile* paramfile);
  void ReInit(const MeshInterface* MP);
};

/*---------------------------------------------------*/

class DwrFem3d : public Q23d
{
 protected:
  typedef Transformation3d<BaseQ13d>         TransQ1;
  FiniteElement<3,2,TransQ1,BaseQ13dPatch>   LowOrderFem;
  HNStructureQ1*                             HNLow;
  
  void TransformationQ1(FemInterface::Matrix& T, int iq) const;

 public:
  DwrFem3d();
  ~DwrFem3d();

  void BasicInit(const ParamFile* paramfile);
  void ReInit(const MeshInterface* MP);
};

/*---------------------------------------------------*/
/*---------------------------------------------------*/

class DwrFemQ1Q22d : virtual public DwrFem2d
{
 protected:
  void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex2d& p0, int i, double s) const;

 public:
  DwrFemQ1Q22d() : DwrFem2d() {}
  ~DwrFemQ1Q22d() {}

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

/*---------------------------------------------------*/

class DwrFemQ1Q23d : virtual public DwrFem3d
{
 protected:
  void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex3d& p0, int i, double s) const;

 public:
  DwrFemQ1Q23d() : DwrFem3d() {}
  ~DwrFemQ1Q23d() {}

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

/*---------------------------------------------------*/
/*---------------------------------------------------*/

class DwrFemQ2Q12d : virtual public DwrFem2d
{
 protected:
  void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex2d& p0, int i, double s) const;

 public:
  DwrFemQ2Q12d() : DwrFem2d() {}
  ~DwrFemQ2Q12d() {}

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

/*---------------------------------------------------*/

class DwrFemQ2Q13d : virtual public DwrFem3d
{
 protected:
  void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex3d& p0, int i, double s) const;

 public:
  DwrFemQ2Q13d() : DwrFem3d() {}
  ~DwrFemQ2Q13d() {}

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
}

/*---------------------------------------------------*/

#endif 
