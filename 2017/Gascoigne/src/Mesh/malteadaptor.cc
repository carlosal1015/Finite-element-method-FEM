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


#include  "malteadaptor.h"
#include  "compareclass.h"
#include  "filescanner.h"
#include  "giota.h"
#include  <fstream>


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
MalteAdaptor::MalteAdaptor(const ParamFile* pf, const DoubleVector& _eta) :
  eta(_eta)
{
  int idim = 0;
  N = eta.size();

  DataFormatHandler DH;

  DH.insert("maxnodes"   ,& maxnodes,10000000);
  DH.insert("dimension"  ,& idim,2);
  DH.insert("coarsening" ,& coarsening,0);
  DH.insert("refining"   ,& refining,1);
  DH.insert("alpha"      ,& alpha,2.);  // lokale Konvergenzordnung bzgl. h
  DH.insert("beta"       ,& beta,1.);   // Skalierung des Aufwandes bzgl. nnodes

  FileScanner FS(DH,pf,"Adaptor");

  if (idim==1) 
  {
    ppp = 1;
  }
  else if (idim==2) 
  {
    ppp = 3;
  }
  else
  {
    ppp=7;
  }

  etasum = accumulate(eta.begin(),eta.end(),0.);
  gamma = pow(0.5,alpha)-1.;

  yfactor = pow(0.5,alpha);
  yfactor *= yfactor;
}

/*-----------------------------------------*/

double MalteAdaptor::Expectation(double theta, double x) const
{
  double neta  = etasum+gamma*theta;
  double nwork = pow(1.+ppp*x,beta);

  return neta * nwork;
}

/*-----------------------------------------*/

double MalteAdaptor::ExpectationCoarsening(double theta, double x) const
{
  double neta  = theta + (etasum-theta) * pow(0.5,alpha);
  double nwork = pow(1.+ppp*x,beta);

  return neta * nwork;
}

/*-----------------------------------------*/

double MalteAdaptor::Expectation(double thetax, double thetay, double x, double y) const
{
  double neta  = thetax*gamma + thetay*(1.-pow(0.5,alpha)) + etasum * pow(0.5,alpha);
  double nwork = x*(1+ppp)+ 1.-y-x + y/(1+ppp);

  nwork = pow(nwork,beta);

  return neta * nwork;
}

/*-----------------------------------------*/

void MalteAdaptor::refine(IntVector& ref) const
{
  if (etasum==0) return;

  int n = eta.size();

  IntVector C(n); 
  iota(C.begin(),C.end(),0);

  typedef CompareObjectBigToSmall<DoubleVector >  CoC;

  sort(C.begin(),C.end(),CoC(eta));
  
  //ofstream file("hyundai");

  double theta  = 0.;
  double minimum = 1.e10;
  double dx     = 1./n;
  int ixmin = 0;
  int ixopt = ixmin;

  int limit = Gascoigne::min_int(n,(maxnodes-N)/(1+ppp));

  for (int i=ixmin; i<limit; i++)
    {
      double x = float(i)*dx;
      x = Gascoigne::min(x,1.);
      theta += eta[C[i]];
      double psi = Expectation(theta,x);
      if (psi<=minimum)
	{
	  minimum = psi;
	  ixopt = i;
	}
    }
  //file.close();

  ref.insert(ref.begin(),C.begin(),C.begin()+ixopt+1);
}

/*-----------------------------------------*/

void MalteAdaptor::coarse(IntVector& coars) const
{
  coars.resize(0);
  if (etasum==0) return;

  int n = eta.size();
  IntVector C(n); 
  iota(C.begin(),C.end(),0);

  typedef CompareObjectBigToSmall<DoubleVector >  CoC;

  sort(C.begin(),C.end(),CoC(eta));
  
  int diff = N-maxnodes;


  cout << "DIFF " << diff << endl;
  if (diff>0)  coars.insert(coars.begin(),C.end()-diff,C.end());
}

/*-----------------------------------------*/

void MalteAdaptor::refine_and_coarse(IntVector& ref, IntVector& coars) const
{
  if (etasum==0) return;

  int n = eta.size();

  IntVector C(n); 
  iota(C.begin(),C.end(),0);

  typedef CompareObjectBigToSmall<DoubleVector >  CoC;

  sort(C.begin(),C.end(),CoC(eta));
  
  double minimum = 1.e10;
  double dx     = 1./n;
  int ixmin = 0;
  int ixopt = ixmin;
  double fac = 0.5;

  DoubleVector  g(n,0.);
  g[0] = eta[C[0]];
  for (int i=1; i<n; i++)
    {
      g[i] = g[i-1] + eta[C[i]];
    }
  int iyopt = 0;
  int zeta = ppp+1;
  
  for (int ix=0; ix<n; ix++)
    {
      double x = float(ix)*dx;
      double y = float(zeta)/(1-zeta) * (maxnodes*fac-n-ix*(zeta-1)) / n;
      
      y = Gascoigne::min(y,1.);
      y = Gascoigne::max(y,0.);
      
      int iy = static_cast<int>(y * n);
      
      double thetax = g[ix];
      double thetay = g[Gascoigne::min_int(n-iy,n-1)];
      double psi = Expectation(thetax,thetay,x,y);
      
      if (psi<=minimum)
	{
	  minimum = psi;
	  ixopt = ix;
	  iyopt = iy;
	}
      if (x+y>1.01) break;
    }
  //cout << "ref,coarse = " << ixopt << " " << iyopt << endl;
  
  coars.insert(coars.begin(),C.end()-iyopt,C.end());

  ref.insert(ref.begin(),C.begin(),C.begin()+ixopt);
}

/*-----------------------------------------*/

void MalteAdaptor::refine(IntVector& ref, IntVector& coars) const
{
  coars.resize(0);
  ref  .resize(0);

  if (!coarsening && !refining) return;
  if (!coarsening && refining)
    {
      refine(ref);
      return;
    }
  if (!refining && coarsening)
    {
      coarse(coars);
      return;
    }
  refine_and_coarse(ref,coars);
}
}

/*-----------------------------------------*/

