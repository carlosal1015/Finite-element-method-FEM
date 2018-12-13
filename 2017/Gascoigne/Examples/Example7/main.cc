/**
*
* Copyright (C) 2008, 2009 by the Gascoigne 3D authors
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


#include  "../Example6/benchmarkproblem.h"
#include  "../Example6/benchmarkmeshagent.h"
#include  "q1lps2d.h"
#include  "multilevelsolver.h"
#include  "dwralgorithm.h"
#include  "numericinterface.h"
#include  "stdsolver.h"
#include  "benchmarkfunctionals.h"
#include  "functionalcontainer.h"
#include  "dirichletdatabycolor.h"

using namespace Gascoigne;

/*----------------------------------------------------------------------------*/

class DualDragProblemDescriptor : public ProblemDescriptorBase
{
 public:
  
  DualDragProblemDescriptor() : ProblemDescriptorBase() {}
  ~DualDragProblemDescriptor() { GetEquationPointer() = NULL; }

  std::string GetName() const { return "DualDragProblemDescriptor";}

  void BasicInit(const ParamFile* filename, const Equation* equation, const ResidualFunctional* fp)
  {
    ProblemDescriptorBase::BasicInit(filename);
    
    GetEquationPointer() = const_cast<Equation*>(equation);

    nvector<int>    comp  = fp->GetComps();
    nvector<double> scale = fp->GetScales();
    std::set<int>   col   = fp->GetColors();
    
    GetDirichletDataPointer()  = new DirichletDataByColor(comp,col,scale);
  }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface*    NewDiscretization(int level) const { return new Q1Lps2d; }
  SolverInterface*            NewSolver(int level)         const { return new StdSolver;}
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
  BenchMarkProblemDescriptor P1;
  P1.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("NavierStokesBenchmark", &P1);

  /////////////
  // Functional
  /////////////
    
  DragFunctional j;   

  ////////////////////
  // Adjoint Problem
  ////////////////////
    
  DualDragProblemDescriptor P2;
  P2.BasicInit(&paramfile,P1.GetEquation(),&j);
  PC.AddProblem("DualDrag", &P2);

  //////////////
  // Discretization etc
  //////////////

  Numeric          N;
  MultiLevelSolver S;
  DwrAlgorithm     B;

  B.BasicInit(&paramfile,&S,&N,&PC);
  B.AdaptiveLoop("NavierStokesBenchmark","DualDrag",j);

  return 0;
}

