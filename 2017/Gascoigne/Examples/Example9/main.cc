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


#include  "meshagent.h"
#include  "q22d.h"
#include  "numericinterface.h"
#include  "problemdescriptorbase.h"
#include  "waveequation.h"
#include  "newmarkalgorithm.h"
#include  "newmarksolver.h"

using namespace Gascoigne;

/* ----------------------------------------- */

class WaveInitialCondition : public DomainInitialCondition
{
 public:
  
  WaveInitialCondition() : DomainInitialCondition() {}
  std::string GetName() const { return "WaveInitialCondition";}
  int        GetNcomp() const { return 1; }

  double operator()(int c, const Vertex2d& v)const 
  {
    //return StiffnessFunction(v.x(),0.5,20.);
    double r = (v.x()-0.5)*(v.x()-0.5)+(v.y()-0.5)*(v.y()-0.5);
    return exp(-50.*r);
  }
};

/*----------------------------------------------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "OceanCrustProblem";}
  void BasicInit(const ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new WaveEquation;
    GetInitialConditionPointer() = new WaveInitialCondition;
    GetBoundaryEquationPointer() = new WaveBoundaryEquation;
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface* NewDiscretization(int level) const { return new Q22d; }
  SolverInterface*         NewSolver(int level)         const { return new NewmarkSolver;}
  MeshAgentInterface*      NewMeshAgent()               const { return new MeshAgent;}
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
  PC.AddProblem("Wave", &LPD);

  //////////////
  // Discretization etc
  //////////////

  Numeric N;

  MultiLevelSolver    S;
  NewmarkAlgorithm    B;
  B.BasicInit(&paramfile,&S,&N,&PC);
  B.Run("Wave");

  return 0;
}

