/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#include  "optadaptor.h"
#include  "compareclass.h"

#include  <fstream>
#include  "giota.h"

using namespace std;

/*********************************************************************/

namespace Gascoigne
{
OptAdaptor::OptAdaptor
(AdaptorData& inf, DoubleVector& e,const DoubleVector& v) 
: info(inf), vol(v), eta(e)
{
  p = info.dim();
  d = 2;
  if(d==2)
    {
      p2 = 4;
      p4 = 16;
    }
  else
    {
      p2 = 8;
      p4 = 64;
    }
  dd = d;
  pp = p;
  factor = 1.;
  
  prepare();
}

/*********************************************************************/

void OptAdaptor::prepare()
{
  n_aimed = Gascoigne::max_int(1,static_cast<int>(info.rfactor()*info.ncells()));
  n_aimed = Gascoigne::min_int(info.maxnodes(),n_aimed);

  info.reset();

  int alpha = static_cast<int>(pp+dd);

  info.eta() = eta.norm_l1();

  double integral = 0.;
  for (int i=0; i<eta.size(); i++)
    {
      if (vol[i]>0)
	{
	  double h  = pow(vol[i],1./dd);
	  eta[i]   /= pow(h,alpha);
	  integral += vol[i] * pow(eta[i],dd/alpha);
	}
    }
  //info.cs() = integral;
  co = pow(integral/n_aimed,1./dd);

  double min = 1.e8, max = 0.;
  for (int i=0; i<eta.size(); i++)
    {
      if (vol[i]>0.)
	{
	  double h       = pow(vol[i],1./dd);
	  double h_aimed = co * pow(eta[i],-1./alpha);
	  double fac = h/h_aimed;
	  eta[i] = fac;
	  min = Gascoigne::min(min,fac);
	  max = Gascoigne::max(max,fac);
	}
    }
  info.minf() = min;
  info.maxf() = max;
  marge = n_aimed - info.ncells();
}

/*********************************************************************/

void OptAdaptor::coarse(IntVector& coarselist)
{
  //  eta = h/h_aimed

  if (info.cfactor()<=0.) return;

  for (int i=0; i<eta.size(); i++)
    {
      if (vol[i]>0.)
	{
	  if (info.cfactor() * eta[i] < 1.)
	    {
	      coarselist.push_back(i);
	    }
	}
    }
  return;

  /*************************** old *******************/

  for (int i=0; i<eta.size(); i++)
    {
      int level = 0; // ???????????
      if ((vol[i]>0) && (level>1))
	{
	  int test = 0;
	  int nchilds = 0;  //  ??????????????
	  for (int ch=0; ch<nchilds; ch++)
	    {
	      int child = 0;//  ??????????????
	      if (info.cfactor() * eta[child] < 1)
		{
		  test++;
		}
	    }
	  if (test==nchilds)
	    {
	      info.nc() += test;
	      marge     += test;
	      for (int ch=0; ch<nchilds; ch++)
		{
		  int child = 0;//  ??????????????
		  coarselist.push_back(child);
		}	      
	    }
	}
    }
}

/*********************************************************************/

void OptAdaptor::refine(IntVector& reflist)
{
  reflist.resize(0);

  IntVector C(eta.size()); iota(C.begin(),C.end(),0);
  sort(C.begin(),C.end(),CompareObjectBigToSmall<DoubleVector > (eta)); 
  
  int i = 0; used = 0;

  while((marge>used) && (i<eta.size()))
    {
      reflist.push_back(C[i]);
      i++; 
      used += 4;
    }

  info.nr() += i;
  marge     -= used;
}

/*********************************************************************/


/*********************************************************************/

void OptAdaptor::RefineGnuplot(IntVector& reflist)
{
  reflist.resize(0);

  IntVector C(eta.size()); 
  iota(C.begin(),C.end(),0);

  typedef CompareObjectBigToSmall<DoubleVector >  CoC;

  sort(C.begin(),C.end(),CoC(eta)); 
  
  int i = 0; used = 0;

  ofstream cmdfile("cmd.gpl");
  cmdfile << "plot  \"eta.dat\" using 1:2 title \"eta\" with lines lw 0" << endl;
  cmdfile << " pause -1" << endl;
  cmdfile << "plot  \"eta.dat\" using 1:3 title \"delta\" with lines lw 0"<<endl;
  cmdfile << " pause -1" << endl;
  cmdfile.close();

  ofstream file("eta.dat");

  for(int ii=0;ii<eta.size();ii++)
    {
      double e = eta[C[ii]];      
      if( (ii>=1) && (e>0) )
	{
 	  double delta = eta[C[ii-1]]-e;
	  file << ii << "\t" << e <<  "\t"<<  delta << endl;
	}
    }

  file.close();
  system("gnuplot cmd.gpl");

  while((marge>used) && (i<eta.size()))
    {
     reflist.push_back(C[i]);
      i++; 
      used += 4;
    }

  info.nr() += i;
  marge     -= used;
}
}

/*********************************************************************/

