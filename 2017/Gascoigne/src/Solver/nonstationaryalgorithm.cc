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


#include  "nonstationaryalgorithm.h"
#include  "stdtimesolver.h"
#include  "compose_name.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
/*-----------------------------------------*/

void NonstationaryAlgorithm::BasicInit(const ParamFile* paramfile, MultiLevelSolver* MLS,
				       const NumericInterface* NI,
				       const ProblemContainer* PC)
{
  MultiLevelAlgorithm::BasicInit(paramfile,MLS,NI,PC);

  DataFormatHandler DFH;
  DFH.insert("dt"    ,&dt    ,1.);
  DFH.insert("theta" ,&theta ,1.);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(paramfile,"Loop");

  time = 0.;
}

/*-------------------------------------------------*/

void NonstationaryAlgorithm::TimeInfoBroadcast()
{
  for (int l=0; l<GetMultiLevelSolver()->nlevels(); l++)
    {
      GetSolver(l)->SetTimeData(dt,theta,time);
    }
}

/*-------------------------------------------------*/

void NonstationaryAlgorithm::InitSolution(const string& initial, VectorInterface& u) const
{
  GetMultiLevelSolver()->GetSolver()->Zero(u);

  if      (initial=="analytic") GetSolver()->SolutionInit(u);
  else if (initial=="boundary") GetSolver()->BoundaryInit(u);
  else if (initial=="zero")     GetSolver()->Zero(u);
  else                          GetSolver()->Read(u,initial);

  GetSolver()->SetBoundaryVector(u);
  GetSolver()->SubtractMean(u);
  GetSolver()->Visu("Results/solve",u,0);
}

/*-----------------------------------------*/

void NonstationaryAlgorithm::ImplicitEuler(const std::string& problemlabel)
{
  theta = 1.;
  ThetaScheme(problemlabel);
}

/*-----------------------------------------*/

void NonstationaryAlgorithm::CrankNicholson(const std::string& problemlabel)
{
  theta = 0.5;
  ThetaScheme(problemlabel);
}

/*-----------------------------------------*/

void NonstationaryAlgorithm::ThetaScheme(const std::string& problemlabel)
{
  int    niter;
  string initial;

  DataFormatHandler DFH;
  DFH.insert("niter", &niter, 1);
  DFH.insert("initial", &initial,"boundary");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");
  
  assert(theta>0.);

  GetMultiLevelSolver()->ReInit(problemlabel);
  GetSolver()->OutputSettings();

  TimeInfoBroadcast();

  VectorInterface u("u"), f("f");
  ReInitVector(u);
  ReInitVector(f);
  InitSolution(initial,u);

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();

  for (int iter=1; iter<=niter; iter++)
    {
      //
      // rhs fuer alten Zeitschritt
      //
      GetSolver()->Zero(f);
      if (theta!=1.) 
	{
	  GetSolver()->Rhs(f,1./theta-1.);
	  const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver());
	  S->StdSolver::Form(f,u,1.-1./theta);
	}
      GetSolver()->MassMatrixVector(f,u,1./(theta*dt));
      
      // neuer Zeitschritt
      //
      time += dt;
      cout << "\n============== " << iter << " ==== theta-scheme === ";
      cout << " [t,dt] "<< time << " " << dt << "\n";

      TimeInfoBroadcast();

      GetSolver()->Rhs(f,1.);
      GetSolver()->SetBoundaryVector(f);
      GetSolver()->SetBoundaryVector(u);

      nlinfo.reset();
  
      Newton(u,f,nlinfo);

      GetSolver()->Visu("Results/solve",u,iter);
      string name = "Results/solve";
      compose_name(name,iter);
      GetSolver()->Write(u,name);
    }
  DeleteVector(u);
  DeleteVector(f);
}

/*-----------------------------------------*/

void NonstationaryAlgorithm::FractionalStepThetaScheme(const std::string& problemlabel)
{
  int    niter;
  string initial;

  DataFormatHandler DFH;
  DFH.insert("niter", &niter, 1);
  DFH.insert("initial", &initial,"boundary");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");
  
  assert(theta>0.);

  GetMultiLevelSolver()->ReInit(problemlabel);
  GetSolver()->OutputSettings();

  TimeInfoBroadcast();

  VectorInterface u("u"), f("f");
  ReInitVector(u);
  ReInitVector(f);
  InitSolution(initial,u);

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();

  double gamma = 1.-sqrt(0.5);
  double alpha = 2.-sqrt(2);  // arbitrary in (0.5,1)
  double c     = 0.;

  for (int iter=1; iter<=niter; iter++)
    {
      int step = (iter-1)%3;

      if      (step==0) { theta = gamma*alpha;              c = (1.-alpha)/alpha;}
      else if (step==1) { theta = (1.-2.*gamma)*(1.-alpha); c = alpha/(1.-alpha);}
      else if (step==2) { theta = gamma*alpha;              c = (1.-alpha)/alpha;}
      //
      // rhs fuer alten Zeitschritt
      //
      GetSolver()->Zero(f);
      if (step!=1) GetSolver()->Rhs(f,1./alpha);

      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver());
      S->StdSolver::Form(f,u,-c);
      GetSolver()->MassMatrixVector(f,u,1./(theta*dt));
      
      // neuer Zeitschritt
      //
      double k = 0.;
      if (step==1) k = dt*(1.-2.*gamma);
      else         k = dt*gamma;
      time += k;

      TimeInfoBroadcast();

      cout << "\n============== " << iter << " ==== fractional-theta === ";
      cout << " [t,dt] "<< time << " " << k << "\n";

      if (step==1) GetSolver()->Rhs(f,1./(1.-alpha));
      GetSolver()->SetBoundaryVector(f);
      GetSolver()->SetBoundaryVector(u);

      nlinfo.reset();
  
      Newton(u,f,nlinfo);

      GetSolver()->Visu("Results/solve",u,iter);
    }
  DeleteVector(u);
  DeleteVector(f);
}

/*-----------------------------------------*/

}
