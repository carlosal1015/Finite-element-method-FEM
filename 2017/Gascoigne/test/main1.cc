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
#include  "equation.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidebyequation.h"
#include  "meshagent.h"
#include  "problemdescriptorbase.h"
#include  "problemcontainer.h"
#include  "functionalcontainer.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */
class LocalEquation : public Equation
{
public:
  string GetName() const {return "Local";}
  int GetNcomp() const {return 1;}
  void OperatorStrong(DoubleVector& b, const FemFunction& U) const {
    b[0] += U[0].m() - U[0].D();
  }
  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const {
    b[0] += U[0].m()*N.m() + U[0].x()*N.x() + U[0].y()*N.y();
  }
  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const {
    A(0,0) +=  M.m()*N.m() +  M.x()*N.x() + M.y()*N.y();
  }
};

/* ----------------------------------------- */
class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile,
		 const ProblemContainer* PC,
		 const FunctionalContainer* FC=0)
      {
	GetMeshAgentPointer() = new MeshAgent;
	
	if(paramfile==NULL) 
	  {
	    int dim=2;
	    int prerefine=3;
	    string inpname("square.inp");
	    GetMeshAgent()->BasicInit(inpname,dim,prerefine,0);
	  }
	
	StdLoop::BasicInit(paramfile, PC, FC);
      }
};

/*---------------------------------------------------*/
class LocalExactSolution : public ExactSolution
{
public:
  LocalExactSolution() : ExactSolution() {}

  string GetName() const {return "LocalExactSolution";}
  double operator()(int c, const Vertex2d& v)const{return v.x()*v.y();}
//   double operator()(int c, const Vertex2d& v)const{return v.x()*v.y()+11.;}
  int GetNcomp() const { return 1; }
};

/*---------------------------------------------------*/
class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  void BasicInit(const ParamFile* pf) {
    GetEquationPointer()      = new LocalEquation;
    GetExactSolutionPointer() = new LocalExactSolution();
    GetRightHandSidePointer() = new RightHandSideByEquation(GetEquation(), GetExactSolution());
    GetDirichletDataPointer() = new DirichletDataByExactSolution(GetExactSolution());

    ProblemDescriptorBase::BasicInit(pf);
    
    GetBoundaryManager()->AddDirichletData(1,0);
    GetBoundaryManager()->AddDirichletData(2,0);
    GetBoundaryManager()->AddDirichletData(3,0);
    GetBoundaryManager()->AddDirichletData(4,0);
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

  ProblemContainer  PC;
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
