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


#include  "dwralgorithm.h"
#include  "stdsolver.h"
#include  "malteadaptor.h"
#include  "meshagent.h"
#include  "energyestimator.h"
#include  "dwrfem.h"

using namespace Gascoigne;
using namespace std;

/*--------------------------------------------------------*/

DiscretizationInterface* DwrAlgorithm::CreateOtherDiscretization() const
{
  DiscretizationInterface* D;

  int dim = GetSolver()->GetMesh()->dimension();
  if (dim==2) D = new DwrFemQ1Q22d;    
  else        D = new DwrFemQ1Q23d;
  
  D->BasicInit(GetSolver()->GetParamfile());
  return D;
}

/*-------------------------------------------------------*/

void DwrAlgorithm::PrimalResidualsHigher(VectorInterface& f, const VectorInterface& u)
{
  GetSolver()->Zero(f);

  // only necessary if z has additional Dirichlet bc compared to u
  //
  GetSolver()->Rhs (f, -0.5);
  GetSolver()->Form(f,u,0.5);

  DiscretizationInterface*  discretization = GetSolver()->GetDiscretization();
  DiscretizationInterface*  D = CreateOtherDiscretization();

  GetSolver()->SetDiscretization(*D,true);
      
  GetSolver()->Rhs (f,   0.5);
  GetSolver()->Form(f,u,-0.5);
//  GetSolver()->SetBoundaryVectorZero(f);
  
  GetSolver()->SetDiscretization(*discretization);
  delete D;
}

/*--------------------------------------------------------*/

void DwrAlgorithm::DualResidualsHigher(VectorInterface& f, const VectorInterface& u, const VectorInterface& z)
{
  GetSolver()->Zero(f);
  // dual problem
  GetSolver()->AddNodeVector("u",u);

  // standard residual
  //
  GetSolver()->Rhs        (f, -0.5);
  GetSolver()->AdjointForm(f,z,0.5);
  //    S.SetBoundaryVectorZero(f);
  GetSolver()->HNDistribute(f);

  // residual respect Q2 test functions
  //
  DiscretizationInterface*  discretization = GetSolver()->GetDiscretization();
  
  DiscretizationInterface* D = CreateOtherDiscretization();
  GetSolver()->SetDiscretization(*D,true);
  
  GetSolver()->Rhs        (f,   0.5);
  GetSolver()->AdjointForm(f,z,-0.5);
  //    S.SetBoundaryVectorZero(f);
  GetSolver()->HNDistribute(f);
  
  GetSolver()->SetDiscretization(*discretization);
  delete D;

  GetSolver()->DeleteNodeVector("u");
}

/*--------------------------------------------------------*/

void DwrAlgorithm::AdaptiveLoop(const std::string& primallabel, 
				const std::string& duallabel, Functional& J) 
{
  int niter;

  DataFormatHandler DFH;
  DFH.insert("niter", &niter, 1);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");

  GetMultiLevelSolver()->ReInit(primallabel);
  GetSolver()->OutputSettings();

  VectorInterface u("u"), f("f"), z("z");
  GlobalVector    uold;

  for (int iter=1; iter<=niter; iter++)
    {
      cout << "\n=========== DwrLoop ============= " << iter << " ======";
      PrintMeshInformation();

      if (iter>1) GetMultiLevelSolver()->ReInit(primallabel);

      ///////////////////////////////
      // Resize vectors
      ///////////////////////////////

      ReInitVector(u);
      ReInitVector(f);
      ReInitVector(z);

      ///////////////////////////////
      // Inititialize primal solution
      ///////////////////////////////

      if (iter==1) GetSolver()->SolutionInit(u);
      else         GetSolver()->InterpolateSolution(u,uold);

      ///////////////////////////////
      // Solve primal problem
      ///////////////////////////////

      NonLinear(u,f,primallabel,iter);
      
      ///////////////////////////////
      // Compute functionals
      ///////////////////////////////

      double j = GetSolver()->ComputeFunctional(f,u,&J);
      cout.precision(8);
      cout << "Functional: " << j << " "; cout << endl;

      ///////////////////////////////
      // Prepare adjoint problem
      ///////////////////////////////

      GetMultiLevelSolver()->AssembleDualMatrix(u);
      GetMultiLevelSolver()->ComputeIlu(u);
      GetMultiLevelSolver()->SetProblem(duallabel);

      GetSolver()->Zero(f);
      GetSolver()->Rhs(f);
      GetSolver()->SubtractMeanAlgebraic(f);
      GetSolver()->SetBoundaryVector(f);
      GetSolver()->SetBoundaryVector(z);

      ///////////////////////////////
      // Solve adjoint problem
      ///////////////////////////////

      CGInfo info(1.e-6,1.e-14,1,100,"DualIt: ");
      LinearSolve(z,f,info);

      GetSolver()->Visu("Results/dual",z,iter);

      ///////////////////////////////
      // DWR Estimator
      ///////////////////////////////

      int n = GetSolver()->GetMesh()->nnodes();
      DoubleVector   eta(n,0.);

      DualResidualsHigher(f,u,z);
      double est =  GetSolver()->ScalarProductWithFluctuations(eta,f,u);
  
      GetMultiLevelSolver()->SetProblem(primallabel);
      PrimalResidualsHigher(f,u);
      est +=  GetSolver()->ScalarProductWithFluctuations(eta,f,z);
      
      double err = j-J.ExactValue();
      double eff = est/err;

      cout.precision(4);
      cout << scientific << "eta,error,eff: " << est << ",  " << err;
      cout.precision(2);
      cout << fixed << ", " << eff << endl;

      if (iter==niter) break;

      ///////////////////////////////
      // Adaptor
      ///////////////////////////////

      IntVector refnodes, coarsenodes;

      MalteAdaptor A(_paramfile,eta);
      A.refine(refnodes,coarsenodes);

      CopyVector(uold,u);

      GetMeshAgent()->refine_nodes(refnodes);
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    }
  DeleteVector(u);
  DeleteVector(f);
}
