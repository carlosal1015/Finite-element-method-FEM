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


#include  "local.h"
#include  "backup.h"
#include  "monitoring.h"
#include  "compose_name.h"

using namespace std;
using namespace Gascoigne;

/*-------------------------------------------------*/

void LocalLoop::run(const std::string& problemlabel)
{
  _iter=1;

  VectorInterface u("u"), f("f"), dat("dat");
  GlobalVector  ualt;

  Monitoring  Moning;

  _clock_newmesh.start();
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(f);
  GetMultiLevelSolver()->ReInitVector(dat,3);

  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  _clock_newmesh.stop();

  GlobalVector& d = GetMultiLevelSolver()->GetSolver()->GetGV(dat);
  string filename("convection.bup");
  ReadBackUpResize(d,filename);

  GetMultiLevelSolver()->AddNodeVector("beta",dat);

  GetMultiLevelSolver()->GetSolver()->OutputSettings();
  InitSolution(u);
  Moning.BasicInit(GetExactValues());

  cout << "\n================== " << _iter << "================";
  cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes();
  cout << " " << GetMeshAgent()->ncells() << endl;
  Moning.SetMeshInformation(_iter,GetMeshAgent()->nnodes(),GetMeshAgent()->ncells());

  Solve(u,f);
  ComputeGlobalErrors(u);

  _clock_functionals.start();
  nvector<double> juh = Functionals(u,f);
  _clock_functionals.stop();

  nvector<double> eta;

  _clock_estimate.start();
  if (_estimator!="none")
  {
    double est = Estimator(eta,u,f);
    Moning.SetSolutionInformation(_iter,juh,est);
  }
  _clock_estimate.stop();
}
