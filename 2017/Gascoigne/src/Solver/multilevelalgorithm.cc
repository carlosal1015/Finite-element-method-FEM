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


#include  "multilevelalgorithm.h"
#include  "stdmultilevelsolver.h"
#include  "energyestimator.h"
#include  "stdsolver.h"
#include  "meshagent.h"
#include  "malteadaptor.h"
#include  "compose_name.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

/*-----------------------------------------*/

void MultiLevelAlgorithm::BasicInit(const ParamFile* paramfile, MultiLevelSolver* MLS, const NumericInterface* NI,
				    const ProblemContainer* PC)
{
  Algorithm::BasicInit(paramfile,NI);

  _S = MLS;

  GetMultiLevelSolver()->BasicInit(NI,GetMeshAgent(),paramfile,PC);

  DataFormatHandler DFH;
  DFH.insert("coarselevel", &_coarselevel, 0);
  DFH.insert("mgomega",     &_mgomega,     1.);
  DFH.insert("mgtype",      &_mgtype,      "V");
  DFH.insert("linearsolve", &_linearsolve, "gmres");
  DFH.insert("gmresmemsize",&_gmresmemsize,100);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(paramfile,"MultiLevelSolver");

  if((_mgtype!="V") && (_mgtype!="W")) abort();
}

/*-------------------------------------------------*/

void MultiLevelAlgorithm::RunLinear(const std::string& problemlabel, const std::string& solver)
{
  VectorInterface u("u"), f("f"), du("du");


  GetMultiLevelSolver()->ReInit(problemlabel);
  if (solver!="unspecified") _linearsolve=solver;
  
  GetSolver()->OutputSettings();
  PrintMeshInformation();

  ReInitVector(u); 
  ReInitVector(f); 
  ReInitVector(du); 

  GetSolver()->Zero(u); 
  GetSolver()->SetBoundaryVector(u);
  
  // Compute RHS
  
  GetSolver()->Zero(f); 
  GetSolver()->Rhs(f); 
  
  // Compute Residual
  
  GetSolver()->HNAverage(u);
  GetSolver()->Form(f,u,-1.);
  GetSolver()->SetBoundaryVectorZero(f);

  // Assemble Matrix and ILU
  
  GetMultiLevelSolver()->AssembleMatrix(u);
  GetMultiLevelSolver()->ComputeIlu(u);
  
  // Solve Linear System

  GetSolver()->Zero(du); 
  
  CGInfo& info = GetSolverInfos()->GetLInfo();

  LinearSolve(du,f,info);

  cout  << endl << "Linear solver " << info.control().status() << endl << endl;
  GetSolver()->Add(u,1.,du);
  GetSolver()->HNZero(u);

  GetSolver()->Visu("Results/multilevel",u,0);

  DeleteVector(u);
  DeleteVector(f);
  DeleteVector(du);
}

/*-----------------------------------------*/

void MultiLevelAlgorithm::AssembleMatrixAndIlu(VectorInterface& u)
{
  GetMultiLevelSolver()->AssembleMatrix(u);
  GetMultiLevelSolver()->ComputeIlu(u);
}

/*-----------------------------------------*/
 
void MultiLevelAlgorithm::VWCycle(vector<double>& res, vector<double>& rw, 
				  int l, int finelevel, int coarselevel, const string& p,
				  VectorInterface& u, VectorInterface& b, VectorInterface& v)
{
  if(l>coarselevel)
    {
      GetSolver(l)->smooth_pre(u,b,v);
      GetSolver(l)->MatrixResidual(v,u,b);
      
      vector<MgInterpolatorInterface*>& I = GetMultiLevelSolver()->GetInterpolatorPointers();

      I[l-1]-> restrict_zero(GetSolver(l-1)->GetGV(b),GetSolver(l)->GetGV(v));
      GetSolver(l-1)->HNDistribute(b);
      GetSolver(l-1)->SetBoundaryVectorZero(b);
      GetSolver(l-1)->Zero(u);
      
      VWCycle(res,rw,l-1,finelevel,coarselevel,p,u,b,v);

      if (p=="W") VWCycle(res,rw,l-1,finelevel,coarselevel,p,u,b,v);

      rw[l] = GetSolver(l-1)->Norm(u);

      GetSolver(l)  -> Zero(v);
      GetSolver(l-1)-> HNAverage(u);
      I[l-1]-> prolongate_add(GetSolver(l)->GetGV(v),GetSolver(l-1)->GetGV(u));
      GetSolver(l-1)-> HNZero(u);
      GetSolver(l)  -> HNZero(v);
      GetSolver(l)  -> SetBoundaryVectorZero(v);
      GetSolver(l)  -> Add(u,_mgomega,v);
      GetSolver(l)  -> smooth_post(u,b,v);
    }
  else
    {
      GetSolver(l)->smooth_exact(u,b,v);
    }
  if ((l>coarselevel) || (l==finelevel))
    {
      GetSolver(l)  -> MatrixResidual(v,u,b);
      res[l] = GetSolver(l)->Norm(v);
    }
}

/*-----------------------------------------*/
 
void MultiLevelAlgorithm::LinearSolve(VectorInterface& du, const VectorInterface& y, 
				      CGInfo& cginfo)
{
  cginfo.reset();

  GetSolver()->HNAverage(du);

  if(_linearsolve == "mg")
  {
    LinearMGSolve(du,y,cginfo);
  }
  else if(_linearsolve == "gmres")
  {
    GmresSolve(du,y,cginfo);
  }
  else
  {
    cerr << "No linear solver for linearsolve \"" << _linearsolve << "\". Use \"mg\" or \"gmres\"." << endl;
    abort();
  }

  GetSolver()->HNZero(du);
  GetSolver()->SubtractMean(du);
}

/*-----------------------------------------*/
 
void MultiLevelAlgorithm::Precondition(VectorInterface& x, VectorInterface& y)
{
  GetSolver()->Equ(x,1.,y);

  CGInfo pinfo;
  pinfo.user().tol()       = 0.;
  pinfo.user().globaltol() = 0.;
  pinfo.user().maxiter()   = 1;
  pinfo.user().printstep() = 0;
  pinfo.user().text()      = "PrecInfo";

  LinearMGSolve(x,y,pinfo);
}

/*-----------------------------------------*/
 
void MultiLevelAlgorithm::LinearMGSolve(VectorInterface& du, const VectorInterface& y, 
					CGInfo& cginfo)
{
  int nl        = GetMultiLevelSolver()->nlevels();
  int finelevel = nl-1;
  int clevel    = min_int(_coarselevel,finelevel);

  for(int level=clevel; level<nl; level++)
    {
      if (GetSolver(level)->DirectSolver()) clevel=level;
    }
  DoubleVector res(nl,0.), rw(nl,0.);

  VectorInterface _mg0("xmg0"), _mg1("xmg1");

  ReInitVector(_mg0); 
  ReInitVector(_mg1); 

  GetSolver()->Equ(_mg0,1.,y);
  
  bool reached = cginfo.check(GetSolver()->Norm(y),0.);

  for(int it=0; !reached; it++)
    {
      string p = _mgtype;
      VWCycle(res,rw,finelevel,finelevel,clevel,p,du,_mg0,_mg1);
      reached = cginfo.check(res[finelevel],rw[finelevel]);
    }
  DeleteVector(_mg0); 
  DeleteVector(_mg1); 
}

/*-----------------------------------------*/

void MultiLevelAlgorithm::NonLinear(VectorInterface& u, VectorInterface& f,
				    const std::string& problemlabel, int iter)
{
  GetSolver()->SubtractMean(u);
  GetSolver()->SetBoundaryVector(u);

  // Compute RHS
  
  GetSolver()->Zero(f); 
  GetSolver()->Rhs(f); 
  GetSolver()->SetBoundaryVector(f);

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();
  nlinfo.reset();
  
  Newton(u,f,nlinfo);

//   GetSolver()->Visu("Results/solve",u,iter);
//   string name = "Results/solve";
//   compose_name(name,iter);
//   GetSolver()->Write(u,name);
}

/*-----------------------------------------*/

void MultiLevelAlgorithm::RunNonLinear(const std::string& problemlabel, int iter)
{
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetSolver()->OutputSettings();
  PrintMeshInformation();
  
  VectorInterface u("u"), f("f");

  ReInitVector(u); 
  ReInitVector(f); 

  GetSolver()->Zero(u); 
  GetSolver()->SolutionInit(u);

  NonLinear(u,f,problemlabel,iter);

  DeleteVector(u);
  DeleteVector(f);
}

/*-----------------------------------------*/

void MultiLevelAlgorithm::GlobalRefineLoop(const std::string& problemlabel)
{
  int niter;

  DataFormatHandler DFH;
  DFH.insert("niter", &niter, 1);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");
  
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetSolver()->OutputSettings();

  VectorInterface u("u"), f("f");
  GlobalVector    uold;

  for (int iter=1; iter<=niter; iter++)
    {
      cout << "\n======================== " << iter << " === GlobalRefineLoop == ";
      cout << GetMeshAgent()->ncells() << " cells" << endl;

      GetMultiLevelSolver()->ReInit(problemlabel);
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

      ReInitVector(u); 
      ReInitVector(f); 

      if (iter==1) GetSolver()->SolutionInit(u);
      else           GetSolver()->InterpolateSolution(u,uold);

      NonLinear(u,f,problemlabel,iter);

      GetSolver()->Visu("Results/solve",u,iter);
      string name = "Results/solve";
      compose_name(name,iter);
      GetSolver()->Write(u,name);

      if (iter<niter) 
	{
	  CopyVector(uold,u);
	  GetMeshAgent()->global_refine(1);
	  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
	}
    }
  DeleteVector(u);
  DeleteVector(f);
}

/*-----------------------------------------*/

void MultiLevelAlgorithm::LocalRefineLoop(const std::string& problemlabel, FunctionalContainer* FC)
{
  int niter;

  DataFormatHandler DFH;
  DFH.insert("niter", &niter, 1);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");
  
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetSolver()->OutputSettings();

  VectorInterface u("u"), f("f");
  GlobalVector    uold;

  for (int iter=1; iter<=niter; iter++)
    {
      cout << "\n======================== " << iter << " === LocalRefineLoop ===" << endl;

      GetMultiLevelSolver()->ReInit(problemlabel);

      ReInitVector(u); 
      ReInitVector(f); 

      if (iter==1) GetSolver()->SolutionInit(u);
      else         GetSolver()->InterpolateSolution(u,uold);
      
      NonLinear(u,f,problemlabel,iter);

      if (FC!=NULL)
	{
	  int n = FC->size();
	  DoubleVector j(n,0.);
	  int i = 0;
	  for (FunctionalContainer::const_iterator it = FC->begin(); it!=FC->end();++it,++i)
	    {
	      j[i] = GetSolver()->ComputeFunctional(f,u,it->second);
	    }
	  cout << "Functionals: ";
	  for (int i=0; i<j.size(); i++) cout << j[i] << " "; cout << endl;
	}
      StdSolver* S = dynamic_cast<StdSolver*>(GetSolver()); assert(S);
      
      MeshAgent* MA = dynamic_cast<MeshAgent*>(GetMeshAgent()); assert(MA);

      S->GetHierarchicalMeshPointer() = MA -> GetHierarchicalMesh();

      EnergyEstimator E(*S);
      DoubleVector    eta;

      double est = E.Estimator(eta,u,f);

      cout << "Estimator: " << est << endl;

      if (iter==niter) break;

      IntVector refnodes, coarsenodes;

      MalteAdaptor A(_paramfile,eta);
      A.refine(refnodes,coarsenodes);

      CopyVector(uold,u);

      GetMeshAgent()->refine_nodes(refnodes);
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    }
}

/*-----------------------------------------*/

}
