/**
*
* Copyright (C) 2008, 2010 by the Gascoigne 3D authors
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


#include  "meshagent.h"
#include  "q1lps2d.h"
#include  "stdsolver.h"
#include  "onelevelalgorithm.h"
#include  "multilevelalgorithm.h"
#include  "numericinterface.h"
#include  "domainmeanfunctional.h"
#include  "problemdescriptorbase.h"
#include  "usefullfunctionsbd.h"
#include  "navierstokeslps2d.h"

using namespace Gascoigne;

/* ----------------------------------------- */

class BenchMarkDirichletData : public DirichletData
{
protected:
  double vmax;
public:
  BenchMarkDirichletData() {
    vmax = 0.3;
  }
  std::string GetName() const {return "Bench";}
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const {

    double y = v.y();

    b.zero();
    if (color!=80)
      {
        double high = 4.1;
        b[1] = vmax * ParabelFunction(y,0.,high);
      }
  }
};

/*---------------------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "NavierStokesBenchmark";}
  void BasicInit(const ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new NavierStokesLps2d(GetParamFile());
    GetDirichletDataPointer() = new BenchMarkDirichletData;
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface*    NewDiscretization(int level) const { return new Q1Lps2d; }
  SolverInterface*            NewSolver(int level)         const { return new StdSolver;}
  MeshAgentInterface*         NewMeshAgent()               const { return new MeshAgent;}
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
  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);
  ProblemContainer PC;
  PC.AddProblem("NavierStokesBenchmark", &LPD);

  //////////////
  // Discretization etc
  //////////////

  MultiLevelSolver    S;
  Numeric             N;

  MultiLevelAlgorithm B;
  B.BasicInit(&paramfile,&S,&N,&PC);
  B.GlobalRefineLoop("NavierStokesBenchmark");

  return 0;
}

