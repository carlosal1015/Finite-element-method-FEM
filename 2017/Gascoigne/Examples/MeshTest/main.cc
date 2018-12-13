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
#include  "q1lps3d.h"
#include  "stdsolver.h"
#include  "stdmultilevelsolver.h"
#include  "onelevelalgorithm.h"
#include  "multilevelalgorithm.h"
#include  "numericinterface.h"
#include  "domainmeanfunctional.h"
#include  "problemdescriptorbase.h"
#include  "usefullfunctionsbd.h"
#include  "navierstokeslps3d.h"
#include  "boundaryfunction.h"
#include  "gascoignemesh3d.h"

using namespace Gascoigne;

/* ----------------------------------------- */

void Output(int iter, const MeshAgent& MA, const StopWatch& stop)
{
  cout << iter << "\t"
       << MA.nnodes() << "\t" << MA.ncells() << "\t"
       << MA.nlevels() << "\t"
       << stop.read() << "\t" << flush;
  system("ps -C MeshTest -o rss | tail -1");
}

/* ----------------------------------------- */

void refine_mesh(MeshAgent& MA,string type)
{
  const GascoigneMesh3d* GM = dynamic_cast<const GascoigneMesh3d*> (MA.GetMesh(0)); 
  
  if (type=="global") MA.global_refine(1);
  else if (type=="local")
    {
      int nn = GM->nnodes();
      double p = 0.045;
      int nr = 1+static_cast<int> (p*nn);
      set<int> rns;
      for (int i=0;i<nr;++i) rns.insert(rand()%nn);

      nvector<int> rn;
      for (set<int>::const_iterator it = rns.begin();it!=rns.end();++it)
	rn.push_back(*it);
      MA.refine_nodes(rn);
    }
}

/* ----------------------------------------- */

int main(int argc, char** argv)
{
  ParamFile paramfile("gascoigne.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  
  cout << "Global Refinement" << endl
       << "=================" << endl;
  
  if (1)
    {
      StopWatch stop;
      stop.start();
      MeshAgent MA;
      MA.BasicInit(&paramfile);
      stop.stop();
      int iter=0;
      Output(iter,MA,stop);
      for (iter=1;iter<6;++iter)
	{
	  stop.start();
	  refine_mesh(MA,"global");
	  stop.stop();
	  Output(iter,MA,stop);
	}
    }


  cout << endl
       << "Local Refinement" << endl
       << "================" << endl;
  for (int K=0;K<3;++K)
    {
      cout << "==== " << K << endl;
      StopWatch stop;
      stop.start();
      MeshAgent MA;
      MA.BasicInit(&paramfile);
      stop.stop();
      int iter=0;
      Output(iter,MA,stop);
      for (iter=1;iter<6;++iter)
	{
	  stop.start();
	  refine_mesh(MA,"local");
	  stop.stop();
	  Output(iter,MA,stop);
	}
    }

  
  
  return 0;
}

