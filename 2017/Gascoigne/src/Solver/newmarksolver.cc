/**
*
* Copyright (C) 2010 by the Gascoigne 3D authors
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


#include "newmarksolver.h"

using namespace Gascoigne;
using namespace std;

/*----------------------------------------------------------------------------*/

void NewmarkSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  FormWithoutMass(gy,gx,d);
  MassMatrixVector(gy,gx,d);
}  

/*----------------------------------------------------------------------------*/

void NewmarkSolver::FormWithoutMass(VectorInterface& gy, const VectorInterface& gx, double d, double s) const
{
  HNAverage(gx);
  HNAverageData();
  
  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  double alpha = 0.25*dt*dt;
  
  GetDiscretization()->Form(GetGV(gy),GetGV(gx),*EQ,d*alpha);
  
  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
    {
      const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
      GetDiscretization()->BoundaryForm(GetGV(gy),GetGV(gx),BM->GetBoundaryEquationColors(),*BE,d*sqrt(alpha)*s);
    }
  HNZero(gx);
  HNZeroData();
  HNDistribute(gy);
  SubtractMeanAlgebraic(gy);
}

/*----------------------------------------------------------------------------*/

void NewmarkSolver::AssembleMatrix(const VectorInterface& gu, double d)
{
  assert(GetMatrix());
  
  double alpha = 0.25*dt*dt;
  
  const GlobalVector& u = GetGV(gu);
  HNAverage(gu);
  HNAverageData();
  
  GetDiscretization()->Matrix(*GetMatrix(),u,*GetProblemDescriptor()->GetEquation(),d*alpha);
  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
    {
      const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
      GetDiscretization()->BoundaryMatrix(*GetMatrix(),u,BM->GetBoundaryEquationColors(),*BE,d*sqrt(alpha));
    }
  DirichletMatrix();
  HNZero(gu);
  HNZeroData();
  GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),_TP,d);
  StdSolver::DirichletMatrix();
}

/*----------------------------------------------------------------------------*/

