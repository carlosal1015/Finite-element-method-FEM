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


#include  "meshagent.h"
#include  "q1lps2d.h"
#include  "stdsolver.h"
#include  "multilevelsolver.h"
#include  "onelevelalgorithm.h"
#include  "multilevelalgorithm.h"
#include  "numericinterface.h"
#include  "domainmeanfunctional.h"
#include  "problemdescriptorbase.h"
#include  "usefullfunctionsbd.h"
#include  "navierstokeslps2d.h"
#include  "boundaryfunction.h"

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

/* ----------------------------------------- */

class RunderKreis : public BoundaryFunction<2>
{
  double   _r;
  Vertex2d _c;
public :

  std::string GetName() const { return "RunderKreis";}
  void BasicInit(Vertex2d c, double r) {
      _c = c;
      _r = r;
    }
  double operator()(const Vertex2d& c) const {
      double r = - _r;
      for (int i=0; i<2; i++)
        {
          double dx = c[i]-_c[i];
          r += dx * dx;
        }
      return r;
    }
};

/*----------------------------------------------------------------------------*/

class CurvedMeshAgent : public MeshAgent
{
 protected:

  RunderKreis RK;

 public:

  CurvedMeshAgent() : MeshAgent()
    {
      double r = 0.25;
      Vertex2d v(2.,2.);
      RK.BasicInit(v,r);

      AddShape(80,&RK);
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
  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);
  ProblemContainer PC;
  PC.AddProblem("NavierStokesBenchmark", &LPD);

  //////////////
  // Discretization etc
  //////////////

  Numeric N;

  MultiLevelSolver    S;
  MultiLevelAlgorithm B;
  B.BasicInit(&paramfile,&S,&N,&PC);
  B.LocalRefineLoop ("NavierStokesBenchmark");

  return 0;
}

