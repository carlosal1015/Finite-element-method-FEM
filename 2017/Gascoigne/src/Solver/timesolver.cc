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


#include "timesolver.h"
#include  "cg.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

void TimeSolver::SetTimeData(double _dt, double _theta, double _time) 
{
  dt    = _dt; 
  theta = _theta;
  time  = _time;
  assert(dt>0.);
  assert(theta>0.);
  GetProblemDescriptor()->SetTime(time,dt);
}

/*-----------------------------------------*/

void TimeSolver::RegisterMatrix()
{
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  int ncomp = EQ->GetNcomp();
  
  if (GetMassMatrixPointer()==NULL)
    GetMassMatrixPointer() = NewMassMatrix(ncomp,_matrixtype);
    
  StdSolver::RegisterMatrix();
}

/*-----------------------------------------*/

void TimeSolver::ReInitMatrix() 
{
  GetDiscretization()->InitFilter(_PF);
  SparseStructure SA;
  GetDiscretization()->Structure(&SA);
  
  GetMatrix()->ReInit(&SA);
  GetIlu()   ->ReInit(&SA);
  
  GetMassMatrix()->ReInit(&SA);
  GetMassMatrix()->zero();
  GetDiscretization()->MassMatrix(*GetMassMatrix()); 
}

/*-----------------------------------------*/

void TimeSolver::AssembleMatrix(const VectorInterface& gu, double d)
{
  StdSolver::AssembleMatrix(gu,d);
  
  double scale = d/(dt*theta);
  GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),_TP,scale);

  StdSolver::DirichletMatrix();
}

/*-----------------------------------------*/
  
void TimeSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  StdSolver::Form(gy,gx,d);
  
  double scale = d/(dt*theta);
  MassMatrixVector(gy,gx,scale);
}

/*-----------------------------------------*/

void TimeSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  ReInitTimePattern(PDX);
  StdSolver::SetProblem(PDX);
}

/*-----------------------------------------*/

void TimeSolver::ReInitTimePattern(const ProblemDescriptorInterface& PDX)
{
  const Equation* EQ = PDX.GetEquation();
  if (EQ) 
    {
      _TP.reservesize(EQ->GetNcomp(),EQ->GetNcomp(),0.);
      EQ->SetTimePattern(_TP);
    }   
}

/*-----------------------------------------*/

void TimeSolver::MassMatrixVector(VectorInterface& gf, const VectorInterface& gu, double d) const
{
        GlobalVector& f = GetGV(gf);
  const GlobalVector& u = GetGV(gu);
  GetMassMatrix()->vmult_time(f,u,_TP,d);
}

/*-----------------------------------------*/

void TimeSolver::InverseMassMatrix(VectorInterface& u, const VectorInterface& f, CGInfo& info)
{
  CG<TimeSolver,VectorInterface> cg(*this);
  cg.solve(u,f,info);
}

/*-----------------------------------------*/

void TimeSolver::precondition(VectorInterface& u, const VectorInterface& f)
{
  Equ(u,1.,f);
  GetMassMatrix()->PrepareJacobi(1.);
  GetMassMatrix()->Jacobi(GetGV(u));
}

/*-----------------------------------------*/

void TimeSolver::cgvmult(VectorInterface& y, const VectorInterface& x, double d) const
{
  MassMatrixVector(y,x,d);
}

/*-----------------------------------------*/

void TimeSolver::L2Projection(VectorInterface& Gu, VectorInterface& Gf)
{
  GlobalVector& u = GetGV(Gu);
  GlobalVector& f = GetGV(Gf);

  TimePattern TP(u.ncomp());
  TP.zero();
  for (int i=0; i<u.ncomp(); i++) TP(i,i) = 1.;

  u.zero();

  bool reached;
  double s = 1.;
  int iter = 0;
  
  GlobalVector g(u.ncomp(),u.n());
  GlobalVector r(u.ncomp(),u.n());
  GlobalVector d(u.ncomp(),u.n());

  assert(u.ncomp()==g.ncomp());
  assert(u.n()==g.n());

  SimpleMatrix *SM = dynamic_cast<SimpleMatrix *>(GetMassMatrix());
  assert(SM);
  SM->PrepareJacobi(s);

  SM->vmult_time(f,u,TP,-s);

  SM->JacobiVector(u);
  SM->Jacobi(f);

  r.equ(1.,f);
  d.equ(1.,f);
  double Res      = r*r;
  double FirstRes = Res;
  //cout << "\t\tpcg " << iter << "\t" << sqrt(Res) << endl;
  
  if (sqrt(Res)<GetSolverData().GetCgMassGlobalTol()) 
  {
    reached = true;
  }
  else
  {
    reached = false;
  }

  while(!reached && iter<GetSolverData().GetCgMassMaxIter())
    {
      iter++;
      g.zero();
      SM->vmult_time_Jacobi(g,d,TP,s);
      double lambda = Res/(g*d);

      u.add(lambda,d);
      r.add(-lambda,g);

      Res = r*r;
      //cout << "\t\tpcg " << iter << "\t" << sqrt(Res) << endl;
      if (Res < GetSolverData().GetCgMassTol() * GetSolverData().GetCgMassTol() * FirstRes || sqrt(Res)<GetSolverData().GetCgMassGlobalTol()) 
      {
        reached = true;
      }
      double betacg = -(r*g)/(d*g);
      d.sequ(betacg,1.,r);
    }

  SM->Jacobi(u);

//   if(iter==GetSolverData().GetCgMassMaxIter())
//   {
//     cout << "too many iterations" << endl;
//   }
//   else
//   {
//     cout << "converged" << endl;
//   }
}

/*-------------------------------------------------------*/


}
