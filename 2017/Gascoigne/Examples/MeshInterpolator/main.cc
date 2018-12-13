/**
*
* Copyright (C) 2005, 2006 by the Gascoigne 3D authors
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


#include  "projectiondescriptor.h"
#include  "stdmultilevelsolver.h"
#include  "meshinterpolator.h"
#include  "solverinfos.h"
#include  "meshagent.h"
#include  "backup.h"
#include  "vtkvisu.h"

using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  string finename = "start2";
  
  ProjectionProblemDescriptor PD;
  PD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("mesh", &PD);
  
  ///////////////////////
  // MESH 
  ///////////////////////

  MeshAgent MA;
  MA.BasicInit(&paramfile);

  ///////////////////////
  // MG Solver
  ///////////////////////

  StdMultiLevelSolver MLS;
  MLS.BasicInit(&MA,&paramfile, &PC);
  MLS.ReInit("mesh");
  MLS.GetSolver()->OutputSettings();
    
  SolverInfos SI;
  SI.BasicInit(&paramfile);
  SI.GetNLInfo().control().matrixmustbebuild() = 1;

  ///////////////////////
  // Rhs Vector
  ///////////////////////

  VectorInterface f("f");
  MLS.ReInitVector(f);
  MLS.Zero(f);
 
  ///////////////////////
  // Mesh Interpolator
  ///////////////////////
    {
      MeshInterpolator MI;
      MI.BasicInit(MLS.GetSolver()->GetDiscretization(),&MA,finename);
      
      GlobalVector uold;
      ReadBackUpResize(uold,finename+".bup");
      MI.RhsForProjection(MLS.GetSolver()->GetGV(f),uold);
    }
  string coarsename = "Results/" + finename + ".projected";

  Monitor moni;
  moni.init(&paramfile,1);
  moni.set_directory("Results");
  MLS.SetMonitorPtr(&moni);

  VectorInterface u("u");
  MLS.ReInitVector(u);
  MLS.Zero(u);

  MLS.Solve(MLS.nlevels()-1,u,f,SI.GetNLInfo());
  MLS.GetSolver()->Visu(coarsename,u,0);

  finename = "start2cell";

  ///////////////////////
  // Mesh Interpolator (CellVector)
  ///////////////////////
    {
      MeshInterpolator MI;
      MI.BasicInit(MLS.GetSolver()->GetDiscretization(),&MA,finename);
      
      GlobalVector uold;
      ReadBackUpResize(uold,finename+".bup");
      MI.InterpolateCellVector(MLS.GetSolver()->GetGV(f),uold);
    }

  coarsename = "Results/" + finename + ".projected";

  VtkVisu V(*MA.GetMesh(0),coarsename,0);
  V.WriteCellData(MLS.GetSolver()->GetGV(f));

  return 0;
}
