/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#include  "stdloop.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainmeanfunctional.h"
#include  "stokesgls2d.h"
#include  "navierstokesgls2d.h"
#include  "zerodirichletdata.h"
#include  "meshagent.h"
#include  "stdmultilevelsolver.h"
#include  "stdsolver.h"
#include  "q1gls2d.h"
#include  "problemdescriptorbase.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */
class LocalEquation : public NavierStokesGls2d
{
 public:
  string GetName() const {return "Local";}
  LocalEquation() : NavierStokesGls2d() {
    _penalty = 0.; _visc = 0.001;
    ST.delta0 = 1.;
    ST.alpha0 = 1.;
    ST.xeta0 = 6.; 
    cerr << "----------------------------------------------\n";
    cerr << _visc << endl;
  }

};

/* ----------------------------------------- */
class LocalBoundaryRightHandSide : public BoundaryRightHandSide
{
 public:
  string GetName() const {return "Local";}
  int GetNcomp() const {return 3;}
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v, const Vertex2d& n, int col) const {

    //b.zero();
    //assert(b.size()==3);

    if (col==1)
      {
	b[1] += -0.005 * n.x() * N.m();
	b[2] += -0.005 * n.y() * N.m();
	return;
      }
    if (col==2)
      {
	b[1] += -0.02 * n.x() * N.m();
	b[2] += -0.02 * n.y() * N.m();
	return;
      }
  }
};

/* ----------------------------------------- */
class LocalSolver : public StdSolver
{
    DiscretizationInterface* NewDiscretization(int dimension, const string& discname) {
    return new Q1Gls2d;
  }
    void BasicInit(int level, const ParamFile* paramfile, const MeshInterface* MP) {
    DoubleVector d(3); d[0]=0.01; d[1]=0.1; d[2]=0.1; 
    _Dat.SetIluModify(d);
    StdSolver::BasicInit(level, paramfile, MP->dimension());
  }  
};

/* ----------------------------------------- */
class LocalMultiLevelSolver : public StdMultiLevelSolver
{
  SolverInterface* NewSolver(int solverlevel) {
    return new LocalSolver;
  }
};
/* ----------------------------------------- */
class LocalLoop : public StdLoop
{
 public:
  void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC,
		 const FunctionalContainer* FC=0)
      {
	GetMeshAgentPointer() = new MeshAgent;
	
	int dim=2;
	int prerefine=0;
	string inpname("mesh3.gup");
	GetMeshAgent()->BasicInit(inpname,dim,prerefine,0);
	
	GetMultiLevelSolverPointer() = new LocalMultiLevelSolver;
	StdLoop::BasicInit(paramfile, PC, FC);
      }
};

/*---------------------------------------------------*/
class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  void BasicInit(const ParamFile* pf) {
    GetEquationPointer()              = new LocalEquation;
    GetDirichletDataPointer()         = new ZeroDirichletData();
    GetBoundaryRightHandSidePointer() = new LocalBoundaryRightHandSide();
    
    ProblemDescriptorBase::BasicInit(pf);
    
    GetBoundaryManager()->AddDirichletData(3,1);
    GetBoundaryManager()->AddDirichletData(3,2);
    GetBoundaryManager()->AddDirichletData(8,1);
    GetBoundaryManager()->AddDirichletData(8,2);
    
    GetBoundaryManager()->AddBoundaryRightHandSide(1);
    GetBoundaryManager()->AddBoundaryRightHandSide(2);
  }
  string GetName() const {return "Local";}
};

/*---------------------------------------------------*/
class LocalDomainFunctional : public virtual AllDomainFunctional
{
 public:
  
  LocalDomainFunctional() : AllDomainFunctional(1,0)
  {
    ExactValue() = 11.25;
  }
  ~LocalDomainFunctional() {}
  
  string GetName() const {
    return "LocalDomain";
  }
};

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  const ParamFile* paramfile(NULL);
  
  if(argc==2)
    {
      paramfile = new ParamFile(argv[1]); 
    }
  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor LPD;
  LPD.BasicInit(paramfile);

  ProblemContainer PC;
  PC.AddProblem("testproblem",&LPD);

  FunctionalContainer FC;
  LocalDomainFunctional j1;
  FC.AddFunctional("domain",&j1);
  
  /////////////
  // Loop
  /////////////
  LocalLoop loop;
  loop.BasicInit(paramfile,&PC,&FC);
  
  /////////////
  // Functionals
  /////////////

    
  loop.run("testproblem");
  
  if(paramfile!=NULL)
    {
      delete paramfile;
    }
  return 0;
}
