/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 by the Gascoigne 3D authors
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


#include  "stdtimesolver.h"
#include  "simplematrix.h"
#include  "compose_name.h"

using namespace std;

/*-------------------------------------------------------------*/
  
namespace Gascoigne
{
StdTimeSolver::StdTimeSolver()
  : StdSolver(), _MMP(NULL), _dt(0.), _theta(0.), _time(0.), _rhs(0.)
{}

/*-------------------------------------------------------------*/
  
StdTimeSolver::~StdTimeSolver()
{
  if (_MMP) { delete _MMP; _MMP=NULL;}
}

/*-------------------------------------------------------*/

string StdTimeSolver::GetName() const
{
  return "StdTimeSolver";
}

/*-------------------------------------------------------------*/
  
void StdTimeSolver::SetTimeData(double dt, double theta, double time, double oldrhs, double newrhs) 
{
  _dt     = dt;
  _theta  = theta;
  _time   = time;
  _rhs[0] = oldrhs;
  _rhs[1] = newrhs;
  if(oldrhs<0)
  {
    _rhs[0] = (1.-theta)/theta;
    _rhs[1] = 1.;
  }

  GetProblemDescriptor()->SetTime(_time,_dt);
}

/*-------------------------------------------------------------*/

GascoigneVisualization* StdTimeSolver::NewGascoigneVisualization() const {
  GascoigneVisualization* p_gv = new GascoigneVisualization;
  p_gv->set_time(_time);
  return p_gv;
}

/*-------------------------------------------------------------*/

void StdTimeSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  const Equation* EQ = PDX.GetEquation();

  if (EQ) 
    {
      GetTimePattern().reservesize(EQ->GetNcomp(),EQ->GetNcomp(),0.);
      EQ->SetTimePattern(GetTimePattern());
    }
  
  StdSolver::SetProblem(PDX);
}

/*-------------------------------------------------------*/

void StdTimeSolver::RegisterMatrix()
{
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  int ncomp = EQ->GetNcomp();

  if (GetMassMatrixPointer()==NULL)
    GetMassMatrixPointer() = NewMassMatrix(ncomp,_matrixtype);
  
  StdSolver::RegisterMatrix();
}

/*-------------------------------------------------------*/

void StdTimeSolver::ReInitMatrix() 
{
  GetDiscretization()->InitFilter(_PF);
  SparseStructure SA;
  GetDiscretization()->Structure(&SA);

  StdSolver::AddPeriodicNodes(&SA);

  GetMatrix()->ReInit(&SA);
  GetIlu()->ReInit(&SA);

  GetMassMatrix()->ReInit(&SA);
  GetMassMatrix()->zero();
  GetDiscretization()->MassMatrix(*GetMassMatrix()); 
}

/*-------------------------------------------------------*/

MatrixInterface* StdTimeSolver::NewMassMatrix(int ncomp, const string& matrixtype)
{
  return new SimpleMatrix;
}

/*-------------------------------------------------------------*/
  
void Gascoigne::StdTimeSolver::InitialCondition(VectorInterface& gf, double d) const
{
  GlobalVector& f = GetGV(gf);
  HNAverageData();

  const Application* IC  = GetProblemDescriptor()->GetInitialCondition();
  const BoundaryInitialCondition* NIC = GetProblemDescriptor()->GetBoundaryInitialCondition();

  if(IC)
    {
       bool done=false;
       const DomainInitialCondition *DRHS = dynamic_cast<const DomainInitialCondition *>(IC);
       if(DRHS)
       {
         GetDiscretization()->Rhs(f,*DRHS,d);
         done = true;
       }
       const DiracInitialCondition *NDRHS = dynamic_cast<const DiracInitialCondition *>(IC);
       if(NDRHS)
       {
         GetDiscretization()->DiracRhs(f,*NDRHS,d);
         done =true;
       }
       if(!done)
       {
         cerr << "InitialCondition should be either of type DomainRightHandSide or DiracRightHandSide!!!" << endl;
         abort();
       }
    }
  if(NIC)
    {
      assert(NIC->GetNcomp()==f.ncomp());
      const BoundaryManager*  BM   = GetProblemDescriptor()->GetBoundaryManager();
      GetDiscretization()->BoundaryRhs(f,BM->GetBoundaryRightHandSideColors(),*NIC,d);	  
    }
  HNZeroData();
  HNDistribute(gf);
}

/*-------------------------------------------------------*/

void StdTimeSolver::TimeRhsOperator(VectorInterface& gf, const VectorInterface& gu) const 
{
  assert(_theta>0.);
  double d = -(1.-_theta)/_theta;
  StdSolver::Form(gf,gu,d);

  if (_dt>0.)
    {
      GlobalVector& f = GetGV(gf);
      const GlobalVector& u = GetGV(gu);
      
      double d = 1./(_dt*_theta);
      GetMassMatrix()->vmult_time(f,u,GetTimePattern(),d);
    }
}

/*-------------------------------------------------------*/

void StdTimeSolver::MassMatrixVector(VectorInterface& gf, const VectorInterface& gu, double d) const
{
  GlobalVector& f = GetGV(gf);
  const GlobalVector& u = GetGV(gu);
  GetMassMatrix()->vmult_time(f,u,GetTimePattern(),d);
}

/*-------------------------------------------------------*/

void StdTimeSolver::TimeRhs(int k, VectorInterface& gf) const
{
  StdSolver::Rhs(gf,_rhs[k-1]);
}

/*-------------------------------------------------------*/

void StdTimeSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  StdSolver::Form(gy,gx,d);

  if (_dt==0.) return;

  assert(_theta>0.);

  double scale = d/(_dt*_theta);

  const GlobalVector& x = GetGV(gx);
  GlobalVector& y = GetGV(gy);

  assert(y.n()==x.n());
  GetMassMatrix()->vmult_time(y,x,GetTimePattern(),scale);
}

/*-------------------------------------------------------*/

void StdTimeSolver::AssembleMatrix(const VectorInterface& gu, double d)
{
  StdSolver::AssembleMatrix(gu,d);

  if (_dt==0.) return;
  assert(_theta>0.);

  double scale = d/(_dt*_theta);
  GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),GetTimePattern(),scale);

  StdSolver::PeriodicMatrix();
  StdSolver::DirichletMatrix();
}

/*-------------------------------------------------------*/

void StdTimeSolver::L2Projection(VectorInterface& Gu, VectorInterface& Gf)
{
  GlobalVector& u = GetGV(Gu);
  GlobalVector& f = GetGV(Gf);

  TimePattern TP(u.ncomp());
  TP.zero();
  for (int i=0; i<u.ncomp(); i++) TP(i,i) = 1.;

  u.zero();
  f.zero();

  InitialCondition(Gf);

  PrecondCGMass(u,f,TP);
}

/*-------------------------------------------------------*/

string StdTimeSolver::PrecondCGMass(GlobalVector& u, GlobalVector& f, const TimePattern& TP, double s)
{
  bool reached;
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

  r.equ(1,f);
  d.equ(1,f);
  double Res = r*r;
  double FirstRes = Res;
  cout << "\t\tpcg " << iter << "\t" << sqrt(Res) << endl;
  
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
      cout << "\t\tpcg " << iter << "\t" << sqrt(Res) << endl;
      if (Res < GetSolverData().GetCgMassTol() * GetSolverData().GetCgMassTol() * FirstRes || sqrt(Res)<GetSolverData().GetCgMassGlobalTol()) 
      {
        reached = true;
      }
      double betacg = -(r*g)/(d*g);
      d.sequ(betacg,1.,r);
    }

  SM->Jacobi(u);

  if(iter==GetSolverData().GetCgMassMaxIter())
  {
    return "too many iterations";
  }
  else
  {
    return "converged";
  }
}

/*-------------------------------------------------------*/

void StdTimeSolver::SetMassMatrix(Gascoigne::MatrixInterface &MM, bool init)
{
  if(init)
  {
    SparseStructure SA;
    GetDiscretization()->Structure(&SA);

    MM.ReInit(&SA);
    MM.zero();
    GetDiscretization()->MassMatrix(MM);
  }
  GetMassMatrixPointer() = &MM;
}
}
