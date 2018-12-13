/**
*
* Copyright (C) 2008, 2010, 2011 by the Gascoigne 3D authors
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


#include  "algorithm.h"
#include  "meshagent.h"
#include  "backup.h"
#include  "adaptordata.h"
#include  "filescanner.h"
#include  "monitoring.h"
#include  <iomanip>
#include  "gostream.h"
#include  "compose_name.h"
#include  "givensrotation.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  Algorithm::Algorithm() : _MA(NULL), _SI(NULL), _paramfile(NULL), _gmresmemsize(100)
{
}

/*-----------------------------------------*/

Algorithm::~Algorithm()
{
  if(_MA  != NULL) {delete _MA; _MA=NULL;}
  if(_SI  != NULL) {delete _SI; _SI=NULL;}
}

/*-----------------------------------------*/

void Algorithm::BasicInit(const ParamFile* paramfile, const NumericInterface* NI)
{
  _paramfile = paramfile;
  _NI = NI;

  GetMeshAgentPointer() = GetNumeric()->NewMeshAgent();
  GetMeshAgent()->BasicInit(_paramfile);

  GetSolverInfosPointer() = new SolverInfos;
  GetSolverInfos()->BasicInit(_paramfile);

  string command("mkdir -p Results");
  system(command.c_str());
}

/*-------------------------------------------------------*/

void Algorithm::PrintMeshInformation() const
{
  cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes();
  cout << " " << GetMeshAgent()->ncells() << endl;
}

/*-----------------------------------------*/

void Algorithm::Newton(VectorInterface& u, const VectorInterface& f, NLInfo& nlinfo)
{
  VectorInterface y("y"), du("du");

  ReInitVector(du); 
  ReInitVector(y); 

  cout << endl << "  ";
  nlinfo.reset();
  CGInfo& cginfo = nlinfo.GetLinearInfo();

  // Newton Residual

  GetSolver()->Equ(y,1.,f);
  GetSolver()->HNAverage(u);
  GetSolver()->Form(y,u,-1.);
  GetSolver()->HNZero(u);
  GetSolver()->SetBoundaryVectorZero(y);
  GetSolver()->SubtractMeanAlgebraic(y);

  double ny = GetSolver()->Norm(y);
  bool reached = nlinfo.check(0,ny,0.);

  for(int it=1; !reached; it++)
    {
      int nm1 = nlinfo.control().newmatrix();
      int nm2 = nlinfo.control().matrixmustbebuild();
      
      if (nm1+nm2!=0)
	{
	  if (nm1 && !nm2) cout << " N";
	  if (!nm1 && nm2) cout << "M ";
          if (nm1 && nm2)  cout << "MN";
	  AssembleMatrixAndIlu(u);
	  nlinfo.control().matrixmustbebuild() = 0;
	}
      else cout << "  ";
      GetSolver()->Zero(du);

      LinearSolve(du,y,cginfo);

      /////////////////////////////////
      //////// Newton update //////////

      double ndu = GetSolver()->Norm(du);

      if ( (cginfo.control().status()=="exploded") || 
	   (ndu>1.e45) || (!(ndu>=0.)) || (ny>1.e45) || (!(ny>=0.)) )
	{
	  nlinfo.control().status()="diverged";
	  cerr << "linear : " << cginfo.control().status() << endl;
	  cerr << "nonlinear : " << ndu << endl;
	}
      else
	{
	  double omega = 0.7;
	  double relax = 1.;
	  
	  for(int iter=0; iter<nlinfo.user().maxrelax(); iter++)
	    {
	      if(iter==0)
		{
		  GetSolver()->Add(u,relax,du);
		}
	      else
		{
		  GetSolver()->Add(u,relax*(omega-1.),du);
		  relax *= omega;
		}
	      GetSolver()->Equ(y,1.,f);
	      GetSolver()->HNAverage(u);
	      GetSolver()->Form(y,u,-1.);
	      GetSolver()->HNZero(u);
	      GetSolver()->SetBoundaryVectorZero(y);
	      GetSolver()->SubtractMeanAlgebraic(y);

	      ny = GetSolver()->Norm(y);

	      string message = nlinfo.check_damping(iter,ny);

	      if (message=="ok")       break;
	      if (message=="continue") continue;
	      if (message=="exploded") 
		{
		  GetSolver()->Add(u,-relax,du);
		  relax = 0.;
		  cout << "Damping exploded !!!!!" << endl;
		  nlinfo.control().status() = "diverged";
		  break;
		}
	    }
	}
      reached = nlinfo.check(it,ny,ndu);
    }
  DeleteVector(y);
  DeleteVector(du);
}

/*-------------------------------------------------------*/

void Algorithm::CopyVector(GlobalVector& dst, VectorInterface& src) const
{
  GetSolver()->HNAverage(src);
  
  int nn = GetSolver()->GetGV(src).n();
  int cc = GetSolver()->GetGV(src).ncomp();

  dst.ncomp() = cc;
  dst.resize(nn);
  
  dst.equ(1.,GetSolver()->GetGV(src));
  
  GetSolver()->HNZero(src);
}

/*-----------------------------------------*/

void Algorithm::GmresSolve(VectorInterface& x, const VectorInterface& b, 
			   CGInfo& info)
{
  // we could implement restarted gmres for maxiter > gmresmemsize
  // for now maxiter = min(gmresmemsize, linear_maxiter)
  int maxiter = Gascoigne::min(_gmresmemsize,info.user().maxiter());

  int minsize = Gascoigne::max(1,Gascoigne::min(5,maxiter));
  vector<VectorInterface> mem;

  int left_precondition = 1;

  for(int i=0; i<minsize; i++)
    {
      std::string s = "gmres";
      compose_name(s,i);
      mem.resize(i+1,s);
      ReInitVector(mem[i]);
    }
  int reached = 0;

  VectorInterface& v = mem[0];
  VectorInterface  p("gmresp");
  ReInitVector(p);
  ReInitVector(v);

  if (left_precondition)
    {
      GetSolver()->residualgmres(p,x,b);
      Precondition(v,p);
    }
  else
    {
      GetSolver()->residualgmres(v,x,b);
    }
  double norm = GetSolver()->Norm(v);

  GetSolver()->Equ(v,1./norm,v);
  GivensRotation   GR(maxiter,norm);

  for (int n = 1; (n<maxiter) && !reached; n++)
    {
      int m = n-1;
      if (n>=mem.size())
        {
	  int i = mem.size();
	  std::string s = "gmres";
	  compose_name(s,i);
	  mem.resize(i+1,s); 
	  ReInitVector(mem[i]);
	}
      VectorInterface& um = mem[m];
      VectorInterface& un = mem[n];
      if (left_precondition)
	{
	  GetSolver()->vmulteq(p,um,1.);
	  Precondition(un,p);
	}
      else
	{
	  GetSolver()->Zero(p);
	  Precondition(p,mem[m]);
	  GetSolver()->vmulteq(un,p,1.);
	}

      for (int i=0 ; i<n ; i++)
	{
	  VectorInterface& ui = mem[i];
 	  double d = GetSolver()->ScalarProduct(un,ui);
	  GR.matrix(i,m) = d;
	  GetSolver()->Add(un,-d,ui);
	}
      double s = GetSolver()->Norm(un);
      GR.matrix(n,m)  = s;  
      GetSolver()->Equ(un,1./s,un);
      double rho = GR.orthogonalization(m);

      reached = info.check(rho,m);
    }
  //
  // Calculate solution
  //
  nvector<double> h = GR.getcoefficients();

  if (left_precondition)
    {
      for (int i=0 ; i<h.size() ; i++)
	GetSolver()->Add(x,h[i],mem[i]);
    }
  else
    {
      GetSolver()->Zero(p);
      for (int i=0 ; i<h.size()  ; i++)
	GetSolver()->Add(p,h[i], mem[i]);
       GetSolver()->Equ(mem[0],0.,mem[0]);
       Precondition(mem[0],p);
       GetSolver()->Add(x,1.,mem[0]);
    }
  //
  // Delete memory
  //
  DeleteVector(p);
  for (int i=0; i<mem.size(); i++)
    {
      DeleteVector(mem[i]);
    }
}

/*-----------------------------------------*/
 
void Algorithm::Precondition(VectorInterface& x, VectorInterface& y)
{
  GetSolver()->Equ(x,1.,y);
}

/*-------------------------------------------------------*/
}
