/**
*
* Copyright (C) 2008 by the Gascoigne 3D authors
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


//#include  "heatproblem.h"
#include  "benchmarkproblem.h"
#include  "benchmarkmeshagent.h"
#include  "q1lps2d.h"
#include  "multilevelsolver.h"
#include  "nonstationaryalgorithm.h"
#include  "numericinterface.h"
#include  "timesolver.h"

using namespace Gascoigne;

/*----------------------------------------------------------------------------*/

class TimeMultiLevelSolver : public MultiLevelSolver
{
public:
  
  SolverInterface* NewSolver(int solverlevel)   {  return new TimeSolver;}
  string GetName() const {  return "TimeMultiLevelSolver";}
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  //DiscretizationInterface*    NewDiscretization(int level) const { return new Q12d; }
  DiscretizationInterface*    NewDiscretization(int level) const { return new Q1Lps2d; }
  SolverInterface*            NewSolver(int level)         const { return new TimeSolver;}
  MeshAgentInterface*         NewMeshAgent()               const { return new CurvedMeshAgent;}
};

/*----------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("gascoigne.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  
  /////////////
  // Equation
  /////////////
  BenchMarkProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);
  ProblemContainer PC;
  PC.AddProblem("NavierStokesBenchmark", &LPD);

  //////////////
  // Discretization etc
  //////////////

  Numeric                N;
  MultiLevelSolver       S;
  NonstationaryAlgorithm B;

  B.BasicInit(&paramfile,&S,&N,&PC);
  //B.ImplicitEuler("NavierStokesBenchmark");
  B.ThetaScheme("NavierStokesBenchmark");
  //B.FractionalStepThetaScheme("NavierStokesBenchmark");

  return 0;
}

