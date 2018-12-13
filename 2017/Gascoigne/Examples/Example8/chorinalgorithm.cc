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


#include  "chorinalgorithm.h"
#include  "stdtimesolver.h"
#include  "compose_name.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

/*-------------------------------------------------*/

void ChorinAlgorithm::Run(const std::string& problem1, const std::string& problem2)
{
  DataFormatHandler DFH;
  DFH.insert("scheme",  &scheme,"Chorin");
  DFH.insert("initial", &initial,"boundary");
  DFH.insert("dt"     , &dt    ,1.);
  DFH.insert("niter"  , &niter ,1);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");

  vproblem = problem1;
  pproblem = problem2;

  GetMultiLevelSolver()->ReInit(vproblem);
  GetMultiLevelSolver()->ReInit(pproblem);

  GetSolver()->OutputSettings();

  TimeInfoBroadcast();

  if      ( scheme == "Chorin" )        { theta=1. ; Chorin()     ; }
  else if ( scheme == "ChorinUzawa" )   { theta=1. ; ChorinUzawa(); }
  else if ( scheme == "VanKan")         { theta=0.5; VanKan()     ; }
  else assert(0);
}

/*-------------------------------------------------*/

void ChorinAlgorithm::VelocityPredictor(VectorInterface& v, VectorInterface& fv, NLInfo& nlinfo, int iter)
{
  GetMultiLevelSolver()->SetProblem(vproblem);
  GetSolver()->SetBoundaryVector(v);
  GetSolver()->SetBoundaryVector(fv);
  nlinfo.reset();
  Newton(v,fv,nlinfo);
  GetSolver()->Visu("Results/predictor",v,iter);
}

/*-------------------------------------------------*/

void ChorinAlgorithm::PressurePoissonProblem(VectorInterface& q, VectorInterface& fp, CGInfo& cginfo)
{
  GetMultiLevelSolver()->SetProblem(pproblem);
  GetSolver()->Zero(fp);
  GetSolver()->Zero(q);
  GetSolver()->Rhs(fp,1./theta);
  LinearSolve(q,fp,cginfo);
}

/*-------------------------------------------------*/

void ChorinAlgorithm::VelocityProjection(VectorInterface& v, VectorInterface& q, 
					 VectorInterface& fv, CGInfo& cginfo, int iter)
{ 
  GetMultiLevelSolver()->SetProblem(vproblem);
  GetSolver()->Zero(fv);
  GetSolver()->AddNodeVector("pressure",q);
  GetSolver()->Rhs(fv,theta*dt);
  GetSolver()->DeleteNodeVector("pressure");
  GetSolver()->MassMatrixVector(fv,v,1.);
  GetSolver()->Zero(v);
  cginfo.reset();
  GetSplittingSolver()->InverseMassMatrix(v,fv,cginfo); //löst das Gleichungssystem Mv=fv
  GetSolver()->Visu("Results/velocity",v,iter);
}

/*-------------------------------------------------*/

void ChorinAlgorithm::PressureUpdate(VectorInterface& p, VectorInterface& q, int iter)
{
  GetMultiLevelSolver()->SetProblem(pproblem);
  GetSolver()->Equ(p,1.,q);
  GetSolver()->Visu("Results/pressure",p,iter);
}

/*-------------------------------------------------*/

void ChorinAlgorithm::Chorin()
{
  GetMultiLevelSolver()->SetProblem(vproblem);
  VectorInterface v("v"), fv("fv");
  ReInitVector(v);
  ReInitVector(fv);
  InitSolution(initial,v);
  GetMultiLevelSolver()->AssembleMatrix(v);

  GetMultiLevelSolver()->SetProblem(pproblem);
  VectorInterface p("p"), fp("fp"), q("q");
  ReInitVector(p);
  ReInitVector(q);
  ReInitVector(fp);
  GetSolver()->AddNodeVector("velocity",v);
  GetSolver()->Zero(p);
  AssembleMatrixAndIlu(p);

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();
  CGInfo& cginfo = GetSolverInfos()->GetLInfo();

  theta = 1.;
  for (int iter=1; iter<=niter; iter++)
  {
    GetMultiLevelSolver()->SetProblem(vproblem);
    GetSolver()->Zero(fv);
    GetSolver()->MassMatrixVector(fv,v,1./dt);

    // neuer Zeitschritt
    //
    time += dt;
    TimeInfoBroadcast();
    
    cout << "\n============== " << iter << " ==== Chorin === ";
    cout << " [t,dt] "<< time << " " << dt << "\n";

    // velocity predictor
    //
    VelocityPredictor(v,fv,nlinfo,iter);

    // pressure Poisson problem
    //
    PressurePoissonProblem(q,fp,cginfo);
    
    // velocity projection
    //
    VelocityProjection(v,q,fv,cginfo,iter);

    // pressure update
    //
    PressureUpdate(p,q,iter);
  }

  GetSolver()->DeleteNodeVector("velocity");

  DeleteVector(v);
  DeleteVector(fv);

  DeleteVector(p);
  DeleteVector(q);
  DeleteVector(fp);
}

/*-------------------------------------------------*/

void ChorinAlgorithm::ChorinUzawa()
{
  GetMultiLevelSolver()->SetProblem(vproblem);
  VectorInterface v("v"), fv("fv"), u("u");
  ReInitVector(v);
  ReInitVector(u);
  ReInitVector(fv);
  InitSolution(initial,v);
  GetMultiLevelSolver()->AssembleMatrix(v);

  GetMultiLevelSolver()->SetProblem(pproblem);
  VectorInterface p("p"), fp("fp"), pold("pold"), q("q");
  ReInitVector(p);
  ReInitVector(q);
  ReInitVector(pold);
  ReInitVector(fp);
  GetSolver()->AddNodeVector("velocity",v);
  GetSolver()->Zero(p);
  GetSolver()->Zero(pold);
  AssembleMatrixAndIlu(p);

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();
  CGInfo& cginfo = GetSolverInfos()->GetLInfo();

  double alpha = 0.5;
  double mu = 1.;

  for (int iter=1; iter<=niter; iter++)
    {
      GetMultiLevelSolver()->SetProblem(vproblem);
      GetSolver()->Zero(fv);
      GetSolver()->AddNodeVector("pressure",p);
      GetSolver()->Rhs(fv,1.);
      GetSolver()->DeleteNodeVector("pressure");
      GetSolver()->AddNodeVector("pressure",pold);
      GetSolver()->Rhs(fv,-1.);
      GetSolver()->DeleteNodeVector("pressure");
      GetSolver()->MassMatrixVector(fv,v,1./dt);

      // neuer Zeitschritt
      //
      time += dt;
      TimeInfoBroadcast();

      cout << "\n============== " << iter << " ==== Chorin-Uzawa === ";
      cout << " [t,dt] "<< time << " " << dt << "\n";

      // velocity predictor
      //
      VelocityPredictor(v,fv,nlinfo,iter);

      // pressure Poisson problem
      //
      PressurePoissonProblem(q,fp,cginfo);

      // velocity projection
      //
      VelocityProjection(v,q,fv,cginfo,iter);

      // pressure update
      //
      GetMultiLevelSolver()->SetProblem(pproblem);
      GetSolver()->Equ(pold,1.,p);
      GetSolver()->Zero(fp);
      GetSolver()->Rhs(fp,alpha*mu*dt);
      GetSolver()->MassMatrixVector(fp,p,1.);
      cginfo.reset();
      GetSplittingSolver()->InverseMassMatrix(p,fp,cginfo); 

      GetSolver()->Visu("Results/pressure",p,iter);
    }
  GetSolver()->DeleteNodeVector("velocity");
  DeleteVector(v);
  DeleteVector(u);
  DeleteVector(fv);
  DeleteVector(pold);
  DeleteVector(p);
  DeleteVector(q);
  DeleteVector(fp);
}

/*-------------------------------------------------*/

void ChorinAlgorithm::VanKan()
{
  GetMultiLevelSolver()->SetProblem(vproblem);
  VectorInterface v("v"), fv("fv");
  ReInitVector(v);
  ReInitVector(fv);
  InitSolution(initial,v);
  GetMultiLevelSolver()->AssembleMatrix(v);

  GetMultiLevelSolver()->SetProblem(pproblem);
  VectorInterface p("p"), fp("fp"), q("q");
  ReInitVector(p);
  ReInitVector(fp);
  ReInitVector(q);
  GetSolver()->AddNodeVector("velocity",v);
  GetSolver()->Zero(p);
  AssembleMatrixAndIlu(p);

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();
  CGInfo& cginfo = GetSolverInfos()->GetLInfo();

  theta = 1.;  // only during first time step

  for (int iter=1; iter<=niter; iter++)
    {
      GetMultiLevelSolver()->SetProblem(vproblem);
      GetSolver()->Zero(fv);
      if (iter>1)
	{
	  GetSolver()->AddNodeVector("pressure",p);
	  GetSolver()->Rhs(fv,1./theta);
	  GetSolver()->DeleteNodeVector("pressure");
	  const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver());
	  S->StdSolver::Form(fv,v,-(1.-theta)/theta);
	}
      GetSolver()->MassMatrixVector(fv,v,1./(theta*dt));

      // neuer Zeitschritt
      //
      time += dt;
      TimeInfoBroadcast();

      cout << "\n============== " << iter << " ==== VanKan === ";
      cout << " [t,dt] "<< time << " " << dt << "\n";
      
      // velocity predictor
      //
      VelocityPredictor(v,fv,nlinfo,iter);

      // pressure Poisson problem
      //
      PressurePoissonProblem(q,fp,cginfo);

      if (iter>1)
	{
	  GetSolver()->Add(p,(1.-theta)/theta,q);
	  GetSolver()->Visu("Results/pressure",p,iter);
	}

      // velocity projection
      //
      VelocityProjection(v,q,fv,cginfo,iter);

      if (iter==1)
	{
	  PressureUpdate(p,q,iter);
	  theta = 0.5;
	}
    }
  GetSolver()->DeleteNodeVector("velocity");
  DeleteVector(v);
  DeleteVector(fv);
  DeleteVector(p);
  DeleteVector(fp);
}

}
