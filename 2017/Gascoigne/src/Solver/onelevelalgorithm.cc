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

#include  "onelevelalgorithm.h"

#include "cg.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

/*-----------------------------------------*/

  void OneLevelAlgorithm::BasicInit(const ParamFile *paramfile,
                                    const NumericInterface *NI,
                                    const ProblemContainer *PC)
  {
    Algorithm::BasicInit(paramfile, NI);

    _S = GetNumeric()->NewSolver();

    GetSolver()->SetDiscretization(*GetNumeric()->NewDiscretization());
    GetSolver()->BasicInit(paramfile, GetMeshAgent()->GetDimension());

    _PC = PC;
  }

/*-----------------------------------------*/

void OneLevelAlgorithm::Precondition(VectorInterface& x, VectorInterface& y)
{
  /*CGInfo pinfo;
  pinfo.user().tol()       = 1.e-12;
  pinfo.user().globaltol() = 1.e-12;
  pinfo.user().maxiter()   = 1;
  pinfo.user().printstep() = 0;
  pinfo.user().text()      = "PrecInfo";

  JacobiSolver(x,y,pinfo);*/

  GetSolver()->Equ(x,1.,y);


}

/*-----------------------------------------*/

void OneLevelAlgorithm::IluSolver(VectorInterface &du,
                                  const VectorInterface &f,
                                  CGInfo &info)
{
  VectorInterface residual("residual");
  ReInitVector(residual);

  bool ok = info.check(GetSolver()->Norm(f), 0.);
  for (int iter = 0; !ok; iter++) {
    GetSolver()->MatrixResidual(residual, du, f);
    double rnorm = GetSolver()->Norm(residual);

    StdSolver *SS = dynamic_cast<StdSolver *>(GetSolver());
    SS->GetIlu()->solve(GetSolver()->GetGV(residual));
    double cnorm = GetSolver()->Norm(residual);

    GetSolver()->Add(du, 1., residual);
    ok = info.check(rnorm, cnorm);
  }
  DeleteVector(residual);
}

/*-----------------------------------------*/

void OneLevelAlgorithm::JacobiSolver(VectorInterface &du,
                                     const VectorInterface &f,
                                     CGInfo &info)
{
  VectorInterface residual("residual");
  ReInitVector(residual);
  StdSolver* S = dynamic_cast<StdSolver *>(GetSolver());
  
  bool ok = info.check(GetSolver()->Norm(f), 0.);
  
  for (int iter = 0; !ok; iter++) {
    S->MatrixResidual(residual, du, f);
    double rnorm = GetSolver()->Norm(residual);

    S->GetMatrix()->Jacobi(GetSolver()->GetGV(residual));
    S->GetMatrix()->Jacobi(GetSolver()->GetGV(residual));
    
    double cnorm = GetSolver()->Norm(residual);

    S->Add(du, 1., residual);
    ok = info.check(rnorm, cnorm);
  }
  DeleteVector(residual);
}


void OneLevelAlgorithm::CGSolver(VectorInterface& x, const VectorInterface& b, CGInfo& info, bool precondition)
{

  VectorInterface g("g"),h("h"),d("d"),Ad("Ad");

  StdSolver* solver = dynamic_cast<StdSolver*>(GetSolver());
  
  solver->RegisterVector(g);
  solver->ReInitVector  (g);
  solver->RegisterVector(h);
  solver->ReInitVector  (h);
  solver->RegisterVector(d);
  solver->ReInitVector  (d);
  solver->RegisterVector(Ad);
  solver->ReInitVector  (Ad);

  solver->Equ(g,1.,b);	
  solver->vmult(g,x,-1.);


  double res = solver->Norm(g);

  if (res==0.) return;
   
  solver->Equ(d,1.,g);
  
  if (precondition) {
    solver->GetMatrix()->PrepareJacobi(1.);
    Precondition(d,g);
  }

  double gh  = solver->ScalarProduct(g,d);

  solver->Equ(d,-1.,d);

  int reached = 0;

  for (int it=0; !reached; it++)
    {
      solver->Zero(Ad);	
      solver->vmult(Ad,d,1.);

      double alpha = gh / solver->ScalarProduct(d,Ad);

      solver->Add(g,alpha,Ad);
      solver->Add(x,-alpha,d);
      res = solver->Norm(g);

      reached = info.check(res);

      if (reached) break;
      
      solver->Equ(h,1.,g);
      if (precondition) {
        Precondition(h,g);
      }
      double ghold = gh;
      gh   = solver->ScalarProduct(g,h);
      double beta = gh/ghold;
      
      //solver->GVsadd(beta,d,-1.,h);
      solver->Equ(d,beta,d);
      solver->Add(d,-1.,h);
    }

  solver->DeleteVector(g);
  solver->DeleteVector(h);
  solver->DeleteVector(d);
  solver->DeleteVector(Ad);
}

/*-----------------------------------------*/

void OneLevelAlgorithm::RunLinear(const std::string& problemlabel,
                                  const std::string& solver)
{
  GetSolver()->NewMesh(GetMeshAgent()->GetMesh(0));

  GetSolver()->SetProblem(*_PC->GetProblem(problemlabel));

  GetSolver()->OutputSettings();
  PrintMeshInformation();

  VectorInterface u("u"), f("f");

  ReInitVector(u);
  ReInitVector(f);

  GetSolver()->Zero(u);

  // Compute RHS
  GetSolver()->Zero(f);
  GetSolver()->Rhs(f);
  GetSolver()->SetBoundaryVector(f);

  // Assemble Matrix
  GetSolver()->RegisterMatrix();
  GetSolver()->ReInitMatrix();
  GetSolver()->MatrixZero();
  GetSolver()->AssembleMatrix(u, 1.);

  // Solve Linear System

  CGInfo &info = GetSolverInfos()->GetLInfo();

  if(solver == "ilu") {
    GetSolver()->ComputeIlu();
    IluSolver(u, f, info);
  
  } else if (solver == "jacobi") {
    StdSolver *SS = dynamic_cast<StdSolver *>(GetSolver());
    SS->GetMatrix()->PrepareJacobi(1.);
    JacobiSolver(u,f,info);
  
  } else if (solver == "cg") { 
    CGSolver(u,f,info,false);
  
  } else if (solver == "gmres") {
    GetSolver()->ComputeIlu();
    GmresSolve(u,f,info);
   
  } else {
    cerr << "OneLevelAlgorithm::RunLinear(): Unknown solver \"" << solver << "\""
         << endl;
    abort();
  }

  cout  << endl << "Linear solver " << info.control().status() << endl << endl;

  GetSolver()->HNZero(u); /* ?! */
  GetSolver()->Visu("Results/onelevel",u,0);

  DeleteVector(u);
  DeleteVector(f);
}

/*-----------------------------------------*/

void OneLevelAlgorithm::AssembleMatrixAndIlu(VectorInterface& u)
{
  GetSolver()->MatrixZero();
  GetSolver()->AssembleMatrix(u,1.);
  GetSolver()->ComputeIlu(u);
}

/*-----------------------------------------*/

void OneLevelAlgorithm::LinearSolve(VectorInterface& du, const VectorInterface& y, CGInfo& cginfo)
{
  cginfo.reset();
  IluSolver(du,y,cginfo);
}

/*-----------------------------------------*/

void OneLevelAlgorithm::RunNonLinear(const std::string& problemlabel)
{
  GetSolver()->NewMesh(GetMeshAgent()->GetMesh(0));

  GetSolver()->SetProblem(*_PC->GetProblem(problemlabel));

  GetSolver()->OutputSettings();
  PrintMeshInformation();

  VectorInterface u("u"), f("f");

  ReInitVector(u);
  ReInitVector(f);

  GetSolver()->Zero(u);
  GetSolver()->SolutionInit(u);
  GetSolver()->SubtractMean(u);
  GetSolver()->SetBoundaryVector(u);

  // Compute RHS

  GetSolver()->Zero(f);
  GetSolver()->Rhs(f);
  GetSolver()->SetBoundaryVector(f);

  // Memory Matrix

  GetSolver()->RegisterMatrix();
  GetSolver()->ReInitMatrix();
  GetSolver()->MatrixZero();

  NLInfo &nlinfo = GetSolverInfos()->GetNLInfo();

  Newton(u, f, nlinfo);

  GetSolver()->Visu("Results/onelevel", u, 0);

  DeleteVector(u);
  DeleteVector(f);
}

} /* namespace Gascoigne */
