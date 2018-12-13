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
#include  "q12d.h"
#include  "stdsolver.h"
#include  "onelevelalgorithm.h"
#include  "multilevelalgorithm.h"
#include  "numericinterface.h"
#include  "problemdescriptorbase.h"
#include  "laplace2d.h"
#include  "zerodirichletdata.h"
#include  "constantrighthandside.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  
  std::string GetName() const { return "egal";}

  void BasicInit(const ParamFile* pf) 
  {
    GetEquationPointer()      = new Laplace2d;
    GetRightHandSidePointer() = new ConstantRightHandSideData(1,0,123.);
    GetDirichletDataPointer() = new ZeroDirichletData;
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface*    NewDiscretization(int level) const { return new Q12d; }
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
  string solver = "jacobi";
  if(argc>=3)
    {
       solver = argv[2];
    }  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);
  ProblemContainer PC;
  PC.AddProblem("laplace", &LPD);

  //////////////
  // Discretization etc
  //////////////

  Numeric N;

  if (solver=="jacobi")
    {
      cout << "=================================" << endl;
      cout << "Algorithm with one level ilu solver:" << endl;
      cout << "=================================" << endl;
      
      OneLevelAlgorithm A;
      
      A.BasicInit(&paramfile,&N,&PC);
      A.RunLinear("laplace");
    }
  else if (solver=="mg")
    {
      cout << "=================================" << endl;
      cout << "Algorithm with Multilevel solver:" << endl;
      cout << "=================================" << endl;
      
      MultiLevelSolver    S;
      MultiLevelAlgorithm B;
      B.BasicInit(&paramfile,&S,&N,&PC);
      B.RunLinear("laplace");
    }
  return 0;
}

/*----------------------------------------------------------------------------*/
