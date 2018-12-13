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


#include  "patchdiscretization.h"
#include  <fstream>
#include  "sparsestructure.h"
#include  "pressurefilter.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
void PatchDiscretization::Structure(SparseStructureInterface* SI) const
{
  SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
  assert(S);

  S->build_begin(n());
  for(int iq=0;iq<GetPatchMesh()->npatches();iq++)
    {
      nvector<int> indices = GetLocalIndices(iq);
      S->build_add(indices.begin(), indices.end());
    }
  S->build_end();  
}

/* ----------------------------------------- */

void PatchDiscretization::Transformation(FemInterface::Matrix& T, int iq) const
{
  int dim = GetPatchMesh()->dimension();
  int ne = GetPatchMesh()->nodes_per_patch();

  nvector<int> indices = *GetPatchMesh()->IndicesOfPatch(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);
  if(dim==2)
    {
      for(int ii=0;ii<ne;ii++)
	{
	  Vertex2d v = GetPatchMesh()->vertex2d(indices[ii]);
	  T(0,ii) = v.x();               
	  T(1,ii) = v.y();
	}
    }
  else if(dim==3)
    {
      for(int ii=0;ii<ne;ii++)
	{
	  Vertex3d v = GetPatchMesh()->vertex3d(indices[ii]);
	  T(0,ii) = v.x();               
	  T(1,ii) = v.y();
	  T(2,ii) = v.z();
	}
    }
}

/* ----------------------------------------- */

void PatchDiscretization::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      //EQ.cell(GetPatchMesh(),iq,__U,__QN);
      GetIntegrator()->Form(EQ,__F,*GetFem(),__U,__QN,__QC);
      LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void PatchDiscretization::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      //EQ.cell(GetPatchMesh(),iq,__U,__QN);
      GetIntegrator()->AdjointForm(EQ,__F,*GetFem(),__U,__QN,__QC);
      LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void PatchDiscretization::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  BE.SetParameterData(__QP);

  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;

      const IntVector& q = *GetMesh()->PatchOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
      for (int i=0; i<q.size(); i++)
        {
          int ip  = q[i];
          int ile = l[i];

          Transformation(T,ip);
          GetFem()->ReInit(T);

          GlobalToLocal(__U,u,ip);

          GetIntegrator()->BoundaryForm(BE,__F,*GetFem(),__U,ile,col,__QN,__QC);
          LocalToGlobal(f,__F,ip,d);
        }
    }
}

/* ----------------------------------------- */

void PatchDiscretization::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      //EQ.cell(GetPatchMesh(),iq,__U,__QN);
      GetIntegrator()->Matrix(EQ,__E,*GetFem(),__U,__QN,__QC);
      LocalToGlobal(A,__E,iq,d);
    }
  
//   ofstream file("MATRIX");
//   A.Write(file);
}

/* ----------------------------------------- */

void PatchDiscretization::BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  BE.SetParameterData(__QP);

  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;

      const IntVector& q = *GetMesh()->PatchOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
      for (int i=0; i<q.size(); i++)
        {
          int ip  = q[i];
	  int ile = l[i];
          
          Transformation(T,ip);
          GetFem()->ReInit(T);

          GlobalToLocal(__U,u,ip);
          GetIntegrator()->BoundaryMatrix(BE,__E,*GetFem(),__U,ile,col,__QN,__QC);
          LocalToGlobal(A,__E,ip,d);
        }
    }
}

/* ----------------------------------------- */

void PatchDiscretization::MassMatrix(MatrixInterface& A) const
{
  nmatrix<double> T;
  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);
      GetIntegrator()->MassMatrix(__E,*GetFem());
      LocalToGlobal(A,__E,iq,1.);
    }
}

/* ----------------------------------------- */

void Gascoigne::PatchDiscretization::MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const
{
  nmatrix<double> T;

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocal(__U,u,iq);
    GetIntegrator()->MassForm(TP,__F,*GetFem(),__U);
    LocalToGlobal(f,__F,iq,s);
  }
}

/* ----------------------------------------- */

void PatchDiscretization::ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const
{
//   const IntegrationFormulaInterface& IF = ErrorFormula();

  int ncomp = u.ncomp();
  err.ncomp() = ncomp;
  err.reservesize(3);
  err = 0.;

  CompVector<double> lerr(ncomp,3); 

  nmatrix<double> T;

  GlobalToGlobalData();
  ES->SetParameterData(__QP);

  for(int iq=0; iq<GetPatchMesh()->npatches(); iq++)
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

void PatchDiscretization::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  RHS.SetParameterData(__QP);

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocalData(iq);
      GetIntegrator()->Rhs(RHS,__F,*GetFem(),__QN,__QC);
      LocalToGlobal(f,__F,iq,s);
    }
}

/* ----------------------------------------- */

void PatchDiscretization::BoundaryRhs(GlobalVector& f, const IntSet& Colors,  const BoundaryRightHandSide& NRHS, double s) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  NRHS.SetParameterData(__QP);
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;

      const IntVector& q = *GetMesh()->PatchOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
      for (int i=0; i<q.size(); i++)
	{
	  int ip  = q[i];
	  int ile = l[i];

	  Transformation(T,ip);
	  GetFem()->ReInit(T);

	  GlobalToLocalData(ip);
	  GetIntegrator()->BoundaryRhs(NRHS,__F,*GetFem(),ile,col,__QN,__QC);
	  LocalToGlobal(f,__F,ip,s);
	}
    }
}

/* ----------------------------------------- */

double PatchDiscretization::compute_element_mean_matrix(int iq, EntryMatrix& E) const
{
  std::cerr << "\"PatchDiscretization::compute_element_mean_matrix\" not written!" << std::endl;
  abort();
}

/* ----------------------------------------- */

void PatchDiscretization::InitFilter(nvector<double>& F) const
{
  PressureFilter* PF = static_cast<PressureFilter*>(&F);
  assert(PF);

  if (!PF->Active()) return;

  PF->ReInit(GetMesh()->nnodes(),GetMesh()->nhanging());
  nmatrix<double> T;

  int nv = GetPatchMesh()->nodes_per_patch();
  EntryMatrix  E(nv,1);

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      double cellsize = GetIntegrator()->MassMatrix(E,*GetFem());
      PF->AddDomainPiece(cellsize);

      nvector<int> ind = *GetPatchMesh()->IndicesOfPatch(iq);

      for(int j=0; j<ind.size(); j++)
	{
	  int jj = ind[j];
	  for(int i=0; i<ind.size(); i++)
	    {
 	      F[jj] += E(i,j,0,0);
	    }
	}
    }
}

/* ----------------------------------------- */

void PatchDiscretization::DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const
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

void PatchDiscretization::DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex2d& p0,int i,double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  Vertex2d Tranfo_p0;
   
  int iq = GetPatchNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "PatchDiscretization::DiracRhsPoint point not found\n";
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

void PatchDiscretization::DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex3d& p0,int i,double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  Vertex3d Tranfo_p0;
   
  int iq = GetPatchNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "PatchDiscretization::DiracRhsPoint point not found\n";
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

void Gascoigne::PatchDiscretization::GlobalToLocalCell(LocalVector& U, const GlobalVector& u, int iq) const
{
  IntVector cells = GetGascoigneMesh()->GetPatchIndexHandler().GetPatch2Cell(iq);
  U.ReInit(u.ncomp(),cells.size());

  for(int i=0;i<cells.size();i++)
    {
    for(int c=0;c<u.ncomp();++c)
      {
        U(i,c) = u(cells[i],c);
      }
    }
}

/* ----------------------------------------- */

double PatchDiscretization::ComputeBoundaryFunctional(const GlobalVector& u, const IntSet& Colors, const BoundaryFunctional& BF) const 
{
  nmatrix<double> T;

  GlobalToGlobalData();
  BF.SetParameterData(__QP);

  double j=0.;
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;

    const IntVector& q = *GetMesh()->PatchOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
    for (int i=0; i<q.size(); i++)
    {
      int ip  = q[i];
      int ile = l[i];

      Transformation(T,ip);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,ip);
      j += GetIntegrator()->ComputeBoundaryFunctional(BF,*GetFem(),ile,col,__U,__QN,__QC);
    }
  }
  return j;
}

/* ----------------------------------------- */

double PatchDiscretization::ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const 
{
  GlobalToGlobalData();
  F.SetParameterData(__QP);
  
  nmatrix<double> T;
  double j=0.;
  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      j += GetIntegrator()->ComputeDomainFunctional(F,*GetFem(),__U,__QN,__QC);
    }
  return j;
}

/* ----------------------------------------- */

double PatchDiscretization::ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const
{
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

double PatchDiscretization::ComputePointValue(const GlobalVector& u, const Vertex2d& p0,int comp) const
{
  Vertex2d Tranfo_p0;

  int iq = GetPatchNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "PatchDiscretization::ComputePointValue point not found\n";
      abort();
    }

  nmatrix<double> T;
  Transformation(T,iq);
  GetFem()->ReInit(T);

  GlobalToLocal(__U,u,iq);
  
  return GetIntegrator()->ComputePointValue(*GetFem(),Tranfo_p0,__U,comp);
}

/* ----------------------------------------- */

double PatchDiscretization::ComputePointValue(const GlobalVector& u, const Vertex3d& p0,int comp) const
{
  Vertex3d Tranfo_p0;

  int iq = GetPatchNumber(p0,Tranfo_p0);
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

void Gascoigne::PatchDiscretization::EvaluateParameterRightHandSide(GlobalVector& f, const DomainRightHandSide& CF, double d) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  CF.SetParameterData(__QP);

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocalData(iq);
      GetIntegrator()->EvaluateCellRightHandSide(__F,CF,*GetFem(),__QN,__QC);

      f.add(d,__F);
    }
}

/* ----------------------------------------- */

void Gascoigne::PatchDiscretization::EvaluateBoundaryParameterRightHandSide(GlobalVector& f,const IntSet& Colors, const BoundaryRightHandSide& CF, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  CF.SetParameterData(__QP);
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->PatchOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
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

}
