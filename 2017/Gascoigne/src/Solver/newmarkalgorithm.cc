/**
*
* Copyright (C) 2010 by the Gascoigne 3D authors
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


#include "newmarkalgorithm.h"
#include "newmarksolver.h"

using namespace Gascoigne;
using namespace std;

/*----------------------------------------------------------------------------*/

void NewmarkAlgorithm::Run(const std::string& problemlabel)
{
  int    niter;
  string initial;
  
  DataFormatHandler DFH;
  DFH.insert("niter", &niter, 1);
  DFH.insert("initial", &initial,"analytic");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");
  
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetSolver()->OutputSettings();
  
  TimeInfoBroadcast();
  
  VectorInterface u("u"), f("f"), uold("uold");
  ReInitVector(u);
  ReInitVector(uold);
  ReInitVector(f);
  InitSolution(initial,u);
  
  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();
  
  double alpha = 0.25*dt*dt;
  
  for (int iter=1; iter<=niter; iter++)
    {
      //
      // rhs fuer alten Zeitschritt
      //
      GetSolver()->Zero(f);
      NewmarkSolver* NMS0 = dynamic_cast<NewmarkSolver*>(GetSolver());
      if (iter==1)
	{
	  GetSolver()->Rhs(f,alpha);	
	  NMS0->FormWithoutMass(f,u,-1.,-1.);
	  GetSolver()->MassMatrixVector(f,u,1.);
	}
      else
	{
	  GetSolver()->Rhs(f,   alpha); // hier eigentlich zum Zeitpunkt t_n-1	
	  GetSolver()->Rhs(f,2.*alpha);	
	  NMS0->FormWithoutMass(f,uold,-1.,-1.);
	  NMS0->FormWithoutMass(f,u   ,-2.,0.);
	  GetSolver()->MassMatrixVector(f,uold,-1.);
	  GetSolver()->MassMatrixVector(f,u   , 2.);
	}
      GetSolver()->Equ(uold,1.,u);
      
      // neuer Zeitschritt
      //
      time += dt;
      cout << "\n============== " << iter << " ==== Newmark-scheme === ";
      cout << " [t,dt] "<< time << " " << dt << "\n";
      
      TimeInfoBroadcast();
      
      GetSolver()->Rhs(f,alpha);
      GetSolver()->SetBoundaryVector(f);
      GetSolver()->SetBoundaryVector(u);
      
      nlinfo.reset();
      
      Newton(u,f,nlinfo);
      
      GetSolver()->Visu("Results/solve",u,iter);
    }
  GetSolver()->DeleteNodeVector("q");
  DeleteVector(u);
  DeleteVector(f);
}

/*----------------------------------------------------------------------------*/
