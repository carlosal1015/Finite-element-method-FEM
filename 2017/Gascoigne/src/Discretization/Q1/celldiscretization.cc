/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2008, 2011 by the Gascoigne 3D authors
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


#include  "celldiscretization.h"
#include  "sparsestructure.h"
#include  "pressurefilter.h"
#include  "gascoignemesh.h"
#include  "columndiagstencil.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"
#include  "simplematrix.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
void CellDiscretization::Structure(SparseStructureInterface* SI) const
{
  SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
  assert(S);

  S->build_begin(n());
  for(int iq=0;iq<GetMesh()->ncells();iq++)
    {
      IntVector indices = GetLocalIndices(iq);
      S->build_add(indices.begin(), indices.end());
    }
  S->build_end();  
}

/* ----------------------------------------- */

void CellDiscretization::Transformation(FemInterface::Matrix& T, int iq) const
{
  int dim = GetMesh()->dimension();
  int ne = GetMesh()->nodes_per_cell(iq);
	
  IntVector indices = GetMesh()->IndicesOfCell(iq);
  assert(ne==indices.size());
	
  T.memory(dim,ne);
  if(dim==2)
    {
		for(int ii=0;ii<ne;ii++)
			{
			Vertex2d v = GetMesh()->vertex2d(indices[ii]);
			T(0,ii) = v.x();               
			T(1,ii) = v.y();
			}
    }
  else if(dim==3)
    {
		for(int ii=0;ii<ne;ii++)
			{
			Vertex3d v = GetMesh()->vertex3d(indices[ii]);
			T(0,ii) = v.x();               
			T(1,ii) = v.y();
			T(2,ii) = v.z();
			}
    }
}

/* ----------------------------------------- */
/* ----------------------------------------- */
/* ----------------------------------------- */

void CellDiscretization::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  EQ.SetParameterData(__QP);
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      //EQ.cell(GetMesh(),iq,__U,__QN);
      GetIntegrator()->Form(EQ,__F,*GetFem(),__U,__QN,__QC);
      LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void CellDiscretization::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  EQ.SetParameterData(__QP);
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      //EQ.cell(GetMesh(),iq,__U,__QN);
      GetIntegrator()->AdjointForm(EQ,__F,*GetFem(),__U,__QN,__QC);
      LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void CellDiscretization::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  BE.SetParameterData(__QP);
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->CellOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalOnBoundary(col);
      for (int i=0; i<q.size(); i++)
        {
          int iq  = q[i];
          int ile = l[i];

          Transformation(T,iq);
          GetFem()->ReInit(T);

          GlobalToLocal(__U,u,iq);

          GetIntegrator()->BoundaryForm(BE,__F,*GetFem(),__U,ile,col,__QN,__QC);
          LocalToGlobal(f,__F,iq,d);
        }
    }
}

/* ----------------------------------------- */

void CellDiscretization::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  EQ.SetParameterData(__QP);
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      //EQ.cell(GetMesh(),iq,__U,__QN);
      GetIntegrator()->Matrix(EQ,__E,*GetFem(),__U,__QN,__QC);
      LocalToGlobal(A,__E,iq,d);
    }
}

/* ----------------------------------------- */

void CellDiscretization::BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  BE.SetParameterData(__QP);
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->CellOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalOnBoundary(col);
      for (int i=0; i<q.size(); i++)
        {
          int iq  = q[i];
          int ile = l[i];
          
          Transformation(T,iq);
          GetFem()->ReInit(T);

          GlobalToLocal(__U,u,iq);
          GetIntegrator()->BoundaryMatrix(BE,__E,*GetFem(),__U,ile,col,__QN,__QC);
          LocalToGlobal(A,__E,iq,d);
        }
    }
}

/* ----------------------------------------- */

void CellDiscretization::MassMatrix(MatrixInterface& A) const
{
  nmatrix<double> T;
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);
      GetIntegrator()->MassMatrix(__E,*GetFem());
      LocalToGlobal(A,__E,iq,1.);
      //      CellDiscretization::LocalToGlobal(A,__E,iq,1.);
    }
}

/* ------------------------------------------ */
void CellDiscretization::BoundaryMassMatrix(MatrixInterface& A, const IntSet& Colors) const
{
  A.zero();
  nmatrix<double> T;
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->CellOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalOnBoundary(col);
      for (int i=0; i<q.size(); i++)
        {
          int iq  = q[i];
          int ile = l[i];
          
          Transformation(T,iq);
          GetFem()->ReInit(T);

	  GetIntegrator()->BoundaryMassMatrix(__E,*GetFem(),ile);
          LocalToGlobal(A,__E,iq,1.);
        }
    }
  //Diagonaleintraege auf 1 setzen, wenn Eintrag noch null, damit A invertierbar ist.
  const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil*>(A.GetStencil());
  assert(ST);
  SimpleMatrix* SM = dynamic_cast<SimpleMatrix*>(&A);
  assert(SM);
  int n = ST->n();
  int pos;
  for(int i = 0; i < n; i++)
  {
    pos = ST->diag(i);
    if(SM->GetValue(pos) == 0)
    {
      SM->GetValue(pos) = 1;
    }
  }
}

/* ----------------------------------------- */

void Gascoigne::CellDiscretization::MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const
{
  nmatrix<double> T;
 
  for(int iq=0;iq<GetMesh()->ncells();++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocal(__U,u,iq);
    GetIntegrator()->MassForm(TP,__F,*GetFem(),__U);
    LocalToGlobal(f,__F,iq,s);
  }
}

/* ----------------------------------------- */

void CellDiscretization::ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const
{
  int ncomp = u.ncomp();
  err.ncomp() = ncomp;
  err.reservesize(3);
  err = 0.;

  GlobalVector lerr(ncomp,3); 
  lerr.zero();

  nmatrix<double> T;
  
  GlobalToGlobalData();
  ES->SetParameterData(__QP);
  
  for(int iq=0; iq<GetMesh()->ncells(); iq++)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);
      GlobalToLocal(__U,u,iq);
      GetIntegrator()->ErrorsByExactSolution(lerr,*GetFem(),*ES,__U,__QN,__QC);

      for(int c=0;c<ncomp;c++)  
	{
	  err(0,c) += lerr(0,c);
	  err(1,c) += lerr(1,c);
	  err(2,c) = Gascoigne::max(err(2,c),lerr(2,c));
	}
    }
  for(int c=0;c<ncomp;c++)  
    {
      err(0,c) = sqrt(err(0,c));
      err(1,c) = sqrt(err(1,c));
    }
}

/* ----------------------------------------- */

void CellDiscretization::AssembleError(GlobalVector& eta, const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const
{
  int ncomp = u.ncomp();
  err.ncomp() = ncomp;
  err.reservesize(3);
  err = 0.;

  eta.ncomp() = 3;
  eta.reservesize(GetMesh()->ncells());

  GlobalVector lerr(ncomp,3); 
  lerr.zero();

  nmatrix<double> T;
  
  GlobalToGlobalData();
  ES->SetParameterData(__QP);
  
  for(int iq=0; iq<GetMesh()->ncells(); iq++)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);
      GlobalToLocal(__U,u,iq);
      GetIntegrator()->ErrorsByExactSolution(lerr,*GetFem(),*ES,__U,__QN,__QC);

      for(int c=0;c<ncomp;c++)  
	{
	  err(0,c) += lerr(0,c);
	  err(1,c) += lerr(1,c);
	  err(2,c) = Gascoigne::max(err(2,c),lerr(2,c));

	  eta(iq,0) += lerr(0,c);
	  eta(iq,1) += lerr(1,c);
	  eta(iq,2) += lerr(2,c);
	}
      eta(iq,0) = sqrt(eta(iq,0));
      eta(iq,1) = sqrt(eta(iq,1));
      eta(iq,2) = sqrt(eta(iq,2));
    }
  for(int c=0;c<ncomp;c++)  
    {
      err(0,c) = sqrt(err(0,c));
      err(1,c) = sqrt(err(1,c));
    }
}

/* ----------------------------------------- */

void CellDiscretization::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  RHS.SetParameterData(__QP);
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocalData(iq);
      GetIntegrator()->Rhs(RHS,__F,*GetFem(),__QN,__QC);
      LocalToGlobal(f,__F,iq,s);
    }
}

/* ----------------------------------------- */

void CellDiscretization::BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  NRHS.SetParameterData(__QP);
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->CellOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalOnBoundary(col);
      for (int i=0; i<q.size(); i++)
	{
	  int iq  = q[i];
	  int ile = l[i];

	  Transformation(T,iq);
	  GetFem()->ReInit(T);

	  GlobalToLocalData(iq);
	  GetIntegrator()->BoundaryRhs(NRHS,__F,*GetFem(),ile,col,__QN,__QC);
	  LocalToGlobal(f,__F,iq,s);
	}
    }
}

/* ----------------------------------------- */

void CellDiscretization::InitFilter(DoubleVector& F) const
{
  PressureFilter* PF = static_cast<PressureFilter*>(&F);
  assert(PF);

  if (!PF->Active()) return;

  PF->ReInit(GetMesh()->nnodes(),GetMesh()->nhanging());
  nmatrix<double> T;
  for(int iq=0; iq<GetMesh()->ncells(); ++iq)
    {
      int nv = GetMesh()->nodes_per_cell(iq);
      EntryMatrix  E(nv,1);

      Transformation(T,iq);
      GetFem()->ReInit(T);

      double cellsize = GetIntegrator()->MassMatrix(E,*GetFem());
      PF->AddDomainPiece(cellsize);

      IntVector ind = GetMesh()->IndicesOfCell(iq);

      for(int i=0;i<ind.size();i++)
 	{
	  for(int j=0;j<ind.size();j++)
	    {
	      F[ind[j]] += E(i,j,0,0);
	    }
      	}
    }
}

/* ----------------------------------------- */

void CellDiscretization::DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const
{
  int dim = GetMesh()->dimension();
  vector<int> comps = DRHS.GetComps();
  int nn = comps.size();

  vector<double> up(nn,0);
 
  if (dim == 2)
    {
      vector<Vertex2d> v2d = DRHS.GetPoints2d();
      assert(nn==v2d.size());
      
      for(int i=0;i<nn;++i)
	{
	  DiracRhsPoint(f,DRHS,v2d[i],i,s);
	}
    }
  else if (dim == 3)
    {
      vector<Vertex3d> v3d = DRHS.GetPoints3d();
      assert(nn==v3d.size());
      for(int i=0;i<nn;++i)
	{
	  DiracRhsPoint(f,DRHS,v3d[i],i,s);
	}
    }
  else
    {
      cerr << "wrong dim = " << dim << endl;
      abort();
    }
}

/* ----------------------------------------- */

void CellDiscretization::DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex2d& p0,int i,double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  Vertex2d Tranfo_p0;
   
  int iq = GetCellNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "CellDiscretization::DiracRhsPoint point not found\n";
      abort();
    }
  
  nmatrix<double> T;
  Transformation(T,iq);
  GetFem()->ReInit(T);
  
  GlobalToLocalData(iq);
  GlobalToGlobalData();
  DRHS.SetParameterData(__QP);

  GetIntegrator()->DiracRhsPoint(__F,*GetFem(),Tranfo_p0,DRHS,i,__QN,__QC);
  LocalToGlobal(f,__F,iq,s);
}

/* ----------------------------------------- */

void CellDiscretization::DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex3d& p0,int i,double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  Vertex3d Tranfo_p0;
   
  int iq = GetCellNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "CellDiscretization::DiracRhsPoint point not found\n";
      abort();
    }
  
  nmatrix<double> T;
  Transformation(T,iq);
  GetFem()->ReInit(T);
  
  GlobalToLocalData(iq);
  GlobalToGlobalData();
  DRHS.SetParameterData(__QP);

  GetIntegrator()->DiracRhsPoint(__F,*GetFem(),Tranfo_p0,DRHS,i,__QN,__QC);
  LocalToGlobal(f,__F,iq,s);
}

/*-----------------------------------------*/

double CellDiscretization::ComputeBoundaryFunctional(const GlobalVector& u, const IntSet& Colors, const BoundaryFunctional& BF) const 
{
  GlobalToGlobalData();
  BF.SetParameterData(__QP);
  
  nmatrix<double> T;
  double j=0.;
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->CellOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalOnBoundary(col);
      for (int i=0; i<q.size(); i++)
        {
          int iq  = q[i];
          int ile = l[i];

          Transformation(T,iq);
          GetFem()->ReInit(T);

          GlobalToLocal(__U,u,iq);
          j += GetIntegrator()->ComputeBoundaryFunctional(BF,*GetFem(),ile,col,__U,__QN,__QC);
        }
    }

  return j;
}

/* ----------------------------------------- */

double CellDiscretization::ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const
{
  GlobalToGlobalData();
  F.SetParameterData(__QP);
  
  nmatrix<double> T;
  double j=0.;
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      j += GetIntegrator()->ComputeDomainFunctional(F,*GetFem(),__U,__QN,__QC);
    }
  return j;
}

/* ----------------------------------------- */

double CellDiscretization::ComputeErrorDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const
{
  GlobalToGlobalData();
  F.SetParameterData(__QP);

  nmatrix<double> T;
  double j=0.;
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      j += GetIntegrator()->ComputeErrorDomainFunctional(F,*GetFem(),__U,__QN,__QC);
    }
  return j;
}

/* ----------------------------------------- */

double CellDiscretization::ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const
{
  FP.SetParameterData(__QP);

  int dim = GetMesh()->dimension();
  vector<int> comps = FP.GetComps();
  int nn = comps.size();

  vector<double> up(nn,0);
 
  if (dim == 2)
    {
      vector<Vertex2d> v2d = FP.GetPoints2d();
      assert(nn==v2d.size());
      
      for(int i=0;i<nn;++i)
	{
	  up[i] = ComputePointValue(u,v2d[i],comps[i]);
	}
    }
  else if (dim == 3)
    {
      vector<Vertex3d> v3d = FP.GetPoints3d();
      assert(nn==v3d.size());
      for(int i=0;i<nn;++i)
	{
	  up[i] = ComputePointValue(u,v3d[i],comps[i]);
	}
    }
  else
    {
      cout << "wronng dimension: dim = " << dim << endl;
      abort();
    }

  return FP.J(up);
}
/* ----------------------------------------- */

double CellDiscretization::ComputePointValue(const GlobalVector& u, const Vertex2d& p0,int comp) const
{
  Vertex2d Tranfo_p0;

  int iq = GetCellNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "CellDiscretization::ComputePointValue point not found\n";
      abort();
    }

  nmatrix<double> T;
  Transformation(T,iq);
  GetFem()->ReInit(T);

  GlobalToLocal(__U,u,iq);

  return GetIntegrator()->ComputePointValue(*GetFem(),Tranfo_p0,__U,comp);
}

/* ----------------------------------------- */

double CellDiscretization::ComputePointValue(const GlobalVector& u, const Vertex3d& p0,int comp) const
{
  Vertex3d Tranfo_p0;

  int iq = GetCellNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "CellDiscretization::ComputePointValue point not found\n";
      abort();
    }

  nmatrix<double> T;
  Transformation(T,iq);
  GetFem()->ReInit(T);

  GlobalToLocal(__U,u,iq);
  
  return GetIntegrator()->ComputePointValue(*GetFem(),Tranfo_p0,__U,comp);
}

/* ----------------------------------------- */

void CellDiscretization::EvaluateCellRightHandSide(GlobalVector& f, const DomainRightHandSide& CF, double d) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  CF.SetParameterData(__QP);

  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocalData(iq);
      GetIntegrator()->EvaluateCellRightHandSide(__F,CF,*GetFem(),__QN,__QC);

      f.add_node(iq,d,0,__F);
    }
}

/* ----------------------------------------- */

void CellDiscretization::EvaluateBoundaryCellRightHandSide(GlobalVector& f,const IntSet& Colors, const BoundaryRightHandSide& CF, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  CF.SetParameterData(__QP);
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->CellOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalOnBoundary(col);
      for (int i=0; i<q.size(); i++)
        {
          int iq  = q[i];
          int ile = l[i];

          Transformation(T,iq);
          GetFem()->ReInit(T);

          GlobalToLocalData(iq);
          GetIntegrator()->EvaluateBoundaryCellRightHandSide(__F,CF,*GetFem(),ile,col,__QN,__QC);

          f.add_node(iq,d,0,__F);
        }
    }
}

/* ----------------------------------------- */

void CellDiscretization::EvaluateParameterRightHandSide(GlobalVector& f, const DomainRightHandSide& CF, double d) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  CF.SetParameterData(__QP);

  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocalData(iq);
      GetIntegrator()->EvaluateCellRightHandSide(__F,CF,*GetFem(),__QN,__QC);

      f.add(d,__F);
    }
}

/* ----------------------------------------- */ 
  
void Gascoigne::CellDiscretization::EvaluateBoundaryParameterRightHandSide(GlobalVector& f,const IntSet& Colors,
    const BoundaryRightHandSide& CF, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  CF.SetParameterData(__QP);
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->CellOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalOnBoundary(col);
      for (int i=0; i<q.size(); i++)
        {
          int iq  = q[i];
          int ile = l[i];

          Transformation(T,iq);
          GetFem()->ReInit(T);

          GlobalToLocalData(iq);
          GetIntegrator()->EvaluateBoundaryCellRightHandSide(__F,CF,*GetFem(),ile,col,__QN,__QC);

          f.add(d,__F);
        }
    }
}

/* ----------------------------------------- */ 
  
void CellDiscretization::InterpolateDomainFunction(GlobalVector& f, const DomainFunction& DF) const
{
  int dim = GetMesh()->dimension();
  f.zero();

  DoubleVector gf;
  gf.resize(DF.GetNcomp());

  const GlobalData& gnd = GetDataContainer().GetNodeData();
  FemData QH;

  if(dim==2)
  {
    for(int r=0; r<GetMesh()->nnodes(); ++r)
    {
      QH.clear();
      GlobalData::const_iterator p=gnd.begin();
      for(; p!=gnd.end(); p++)
      {
        FemFunction& UH = QH[p->first];
        const GlobalVector& U = *p->second;
        UH.resize(U.ncomp());
        for (int c=0; c<UH.size(); c++)
        {
          UH[c].zero();
          UH[c].m() = U(r,c);
        }
      }

      DF.SetFemData(QH);
      Vertex2d v = GetMesh()->vertex2d(r);
      DF.F(gf,v);
      f.add_node(r,1.,gf);
    }
  }
  else
  {
    for(int r=0; r<GetMesh()->nnodes(); ++r)
    {
      QH.clear();
      GlobalData::const_iterator p=gnd.begin();
      for(; p!=gnd.end(); p++)
      {
        FemFunction& UH = QH[p->first];
        const GlobalVector& U = *p->second;
        for (int c=0; c<UH.size(); c++)
        {
          UH[c].zero();
          UH[c].m() = U(r,c);
        }
      }

      DF.SetFemData(QH);
      Vertex3d v = GetMesh()->vertex3d(r);
      DF.F(gf,v);
      f.add_node(r,1.,gf);
    }
  }
}

/* ----------------------------------------- */  

void CellDiscretization::InterpolateCellDomainFunction(GlobalVector& f, const DomainFunction& DF) const
{
  int dim = GetMesh()->dimension();
  f.zero();

  DoubleVector gf;
  gf.resize(DF.GetNcomp());

  if(dim==2)
  {
    Vertex2d v;
    for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      v.zero();
      for (int in=0; in<4; in++)
      {
	int r = GetMesh()->vertex_of_cell(iq,in);
	v += GetMesh()->vertex2d(r);
      }
      v *= 0.25;
      DF.F(gf,v);
      f.add_node(iq,1.,gf);
    }
  }
  else
  {
    Vertex3d v;
    for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      v.zero();
      for (int in=0; in<8; in++)
      {
	int r = GetMesh()->vertex_of_cell(iq,in);
	v += GetMesh()->vertex3d(r);
      }
      v *= 0.125;
      DF.F(gf,v);
      f.add_node(iq,1.,gf);
    }
  }
}

/* ----------------------------------------- */

void CellDiscretization::Transformation_HM(FemInterface::Matrix& T, const HierarchicalMesh* HM, int iq) const
{
  int dim = GetMesh()->dimension();
  int ne = GetMesh()->nodes_per_cell(iq);

  IntVector indices = HM->GetVertices(iq);
  assert(ne==indices.size());
  swapIndices(indices);

  T.memory(dim,ne);
  if(dim==2)
    {
      for(int ii=0;ii<ne;ii++)
	{
	  Vertex2d v = GetMesh()->vertex2d(indices[ii]);
	  T(0,ii) = v.x();               
	  T(1,ii) = v.y();
	}
    }
  else if(dim==3)
    {
      for(int ii=0;ii<ne;ii++)
	{
	  Vertex3d v = GetMesh()->vertex3d(indices[ii]);
	  T(0,ii) = v.x();               
	  T(1,ii) = v.y();
	  T(2,ii) = v.z();
	}
    }
}

/* ----------------------------------------- */

void CellDiscretization::GlobalToLocal_HM(LocalVector& U, const GlobalVector& u, const HierarchicalMesh* HM, int iq) const
{
  IntVector indices = HM->GetVertices(iq);
  swapIndices(indices);

  U.ReInit(u.ncomp(),indices.size());
  for(int ii=0; ii<indices.size(); ii++) 
    {
      int i = indices[ii];
      U.equ_node(ii,i,u);
    }
}

/* ----------------------------------------- */

void CellDiscretization::swapIndices(IntVector& indices) const
{
  assert(indices.size()>=4);

  int help = indices[2];
  indices[2] = indices[3];
  indices[3] = help;
  if (indices.size()==8)
    {
      help = indices[6];
      indices[6] = indices[7];
      indices[7] = help;
    }
}

/* ----------------------------------------- */

void CellDiscretization::GetVolumes(DoubleVector& a) const
{
  a.resize(GetMesh()->ncells());
  nmatrix<double> T;
  int dim = GetMesh()->dimension();
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);
      if(dim==2)
	{
	  Vertex2d xi;
	  xi.x() = 0.5;
	  xi.y() = 0.5;
	  GetFem()->point(xi);
	}
      else
	{
	  Vertex3d xi;
	  xi.x() = 0.5;
	  xi.y() = 0.5;
	  xi.z() = 0.5;
	  GetFem()->point(xi);
	}
      a[iq] = GetFem()->J();
    }
}
/* ----------------------------------------- */

void CellDiscretization::GetAreas(DoubleVector& a, const IntSet& Colors) const
{
  a.resize(GetMesh()->ncells(),1.);
  nmatrix<double> T;
  int dim = GetMesh()->dimension();
   
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;
    const IntVector& q = *GetMesh()->CellOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalOnBoundary(col);
    for (int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ile = l[i];
      
      Transformation(T,iq);
      GetFem()->ReInit(T);
      if(dim==2)
      {
	Vertex1d xi;
	xi.x() = 0.5;
	GetFem()->point_boundary(ile,xi);
      }
      else
      {
	Vertex2d xi;
	xi.x() = 0.5;
	xi.y() = 0.5;
	GetFem()->point_boundary(ile,xi);
      }
      if(a[iq] == 1.)
	a[iq] = GetFem()->G();
      else
	a[iq] += GetFem()->G();
    }
  }
}
/* ----------------------------------------- */

void CellDiscretization::GetMassDiag(DoubleVector& a) const
{
    a.resize(GetMesh()->nnodes());
    nmatrix<double> T;
    DoubleVector F;
    for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GetIntegrator()->IntegrateMassDiag(F,*GetFem());
	
	//Auf Globalen Vektor verteielen
	IntVector indices = GetLocalIndices(iq);
	for(int ii=0; ii<indices.size(); ii++) 
	{
	    int i = indices[ii];
	    a[i] += F[ii];
	}
    }
}

/* ----------------------------------------- */

void CellDiscretization::GetBoundaryMassDiag(DoubleVector& a) const
{
    a.resize(GetMesh()->nnodes());
    a.equ(1.);
    
    const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
    assert(GMP);
    std::set<int> Colors =  GMP->GetColors();
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
	int col = *p;
	const IntVector& bv = *GMP->VertexOnBoundary(col);
	for(int i=0;i<bv.size();i++)
	{
	    a[bv[i]] = 0;
	}  
    }
    
    DoubleVector F;
    nmatrix<double> T;
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->CellOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalOnBoundary(col);
      for (int i=0; i<q.size(); i++)
        {
          int iq  = q[i];
          int ile = l[i];

          Transformation(T,iq);
          GetFem()->ReInit(T);

	  GetIntegrator()->IntegrateBoundaryMassDiag(F,*GetFem(),ile,col);

          //Auf Globalen Vektor verteielen
	  IntVector indices = GetLocalIndices(iq);
	  for(int ii=0; ii<indices.size(); ii++) 
	  {
	      int i = indices[ii];
	      a[i] += F[ii];
	  }
        }
    }
}

}
