/**
*
* Copyright (C) 2006, 2007, 2008 by the Gascoigne 3D authors
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


#include "q4.h"
#include "pressurefilter.h"
#include "sparsestructure.h"

#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_set>
#define HASHSET std::tr1::unordered_set
#else
#include  <ext/hash_set>
#define HASHSET __gnu_cxx::hash_set
#endif


using namespace std;

namespace Gascoigne
{
/**********************************************************/

void Q4::GlobalToLocalCell(LocalVector& U, const GlobalVector& u, int iq) const
{
  IntVector cells = GetGascoigneMesh()->GetPatchIndexHandler().GetQ4Patch2Cell(iq);
  U.ReInit(u.ncomp(),cells.size());

  for(int i=0;i<cells.size();i++)
  {
    for(int c=0;c<u.ncomp();++c)
    {
      U(i,c) = u(cells[i],c);
    }
  }
}

/**********************************************************/

void Q4::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  HN->CondenseHanging(E,indices);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,E,s);
}

/**********************************************************/

void Q4::Transformation(FemInterface::Matrix& T, int iq) const
{
  int dim = GetPatchMesh()->dimension();
  int ne = GetPatchMesh()->nodes_per_q4patch();

  nvector<int> indices = *GetPatchMesh()->IndicesOfQ4Patch(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);
  if(dim==2)
  {
    for(int ii=0; ii<ne; ii++)
    {
      Vertex2d v = GetPatchMesh()->vertex2d(indices[ii]);
      T(0,ii) = v.x();
      T(1,ii) = v.y();
    }
  }
  else if(dim==3)
  {
    for(int ii=0; ii<ne; ii++)
    {
      Vertex3d v = GetPatchMesh()->vertex3d(indices[ii]);
      T(0,ii) = v.x();
      T(1,ii) = v.y();
      T(2,ii) = v.z();
    }
  }
}

/**********************************************************/

void Q4::ReInit(const MeshInterface* MP)
{
  PatchDiscretization::ReInit(MP);
  HN->ReInit(MP);
}

/**********************************************************/

void Q4::Structure(SparseStructureInterface* SI) const
{
  SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
  assert(S);

  S->build_begin(n());
  for(int iq=0; iq<GetPatchMesh()->nq4patches(); iq++)
  {
    nvector<int> indices = GetLocalIndices(iq);
    HN->CondenseHanging(indices);
    S->build_add(indices.begin(), indices.end());
  }
  HN->SparseStructureDiag(S);
  S->build_end();
}

/**********************************************************/

void Q4::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); iq++)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocal(__U,u,iq);
    //EQ.cell(GetPatchMesh(),iq,__U,__QN);
    GetIntegrator()->Form(EQ,__F,*GetFem(),__U,__QN,__QC);
    PatchDiscretization::LocalToGlobal(f,__F,iq,d);
  }
}

/**********************************************************/

void Q4::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); iq++)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocal(__U,u,iq);
    //EQ.cell(GetPatchMesh(),iq,__U,__QN);
    GetIntegrator()->AdjointForm(EQ,__F,*GetFem(),__U,__QN,__QC);
    PatchDiscretization::LocalToGlobal(f,__F,iq,d);
  }
}

/**********************************************************/

void Q4::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
{
  int dim = GetPatchMesh()->dimension();
  nmatrix<double> T;

  GlobalToGlobalData();
  BE.SetParameterData(__QP);

  /// die cell2patch - liste

  const nvector<IntVector>& q4patch2cell = GetGascoigneMesh()->GetPatchIndexHandler().GetAllQ4Patch2Cell();

  nvector<int> cell2q4patch(GetMesh()->ncells());
  for(int p=0; p<q4patch2cell.size(); p++)
  {
    for(int i=0; i<q4patch2cell[p].size(); i++)
    {
      cell2q4patch[q4patch2cell[p][i]] = p;
    }
  }

  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;

    HASHSET<int> habschon;

    const IntVector& q = *GetMesh()->PatchOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
    for(int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2q4patch[iq];
      int ile = l[i];

      // gabs den patch schon?
      if(habschon.find((ip<<dim)+ile)!=habschon.end())
      {
        continue;
      }
      habschon.insert((ip<<dim)+ile);

      Transformation(T,ip);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,ip);

      GetIntegrator()->BoundaryForm(BE,__F,*GetFem(),__U,ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,d);
    }
  }
}

/**********************************************************/

void Q4::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocal(__U,u,iq);
    //EQ.cell(GetPatchMesh(),iq,__U,__QN);
    GetIntegrator()->Matrix(EQ,__E,*GetFem(),__U,__QN,__QC);
    LocalToGlobal(A,__E,iq,d);
  }

  HN->MatrixDiag(u.ncomp(),A);
//   ofstream file("MATRIX");
//   A.Write(file);
}

/**********************************************************/

void Q4::BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
{
  int dim = GetPatchMesh()->dimension();
  nmatrix<double> T;

  GlobalToGlobalData();
  BE.SetParameterData(__QP);

  /// die cell2patch - liste

  const nvector<IntVector>& q4patch2cell = GetGascoigneMesh()->GetPatchIndexHandler().GetAllQ4Patch2Cell();

  nvector<int> cell2q4patch(GetMesh()->ncells());

  for(int p=0; p<q4patch2cell.size(); p++)
  {
    for(int i=0; i<q4patch2cell[p].size(); i++)
    {
      cell2q4patch[q4patch2cell[p][i]] = p;
    }
  }

  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;
    HASHSET<int> habschon;

    const IntVector& q = *GetMesh()->PatchOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
    for(int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2q4patch[iq];
      int ile = l[i];
	    
      // gabs den patch schon?
      if(habschon.find((ip<<dim)+ile)!=habschon.end())
      {
        continue;
      }
      habschon.insert((ip<<dim)+ile);

      Transformation(T,ip);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,ip);
      GetIntegrator()->BoundaryMatrix(BE,__E,*GetFem(),__U,ile,col,__QN,__QC);
      LocalToGlobal(A,__E,ip,d);
    }
  }
}

/**********************************************************/

void Q4::MassMatrix(MatrixInterface& A) const
{
  nmatrix<double> T;
  for(int iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);
    GetIntegrator()->MassMatrix(__E,*GetFem());
    LocalToGlobal(A,__E,iq,1.);
  }

  HN->MatrixDiag(1,A);
}

/**********************************************************/

void Q4::MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const
{
  nmatrix<double> T;

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocal(__U,u,iq);
    GetIntegrator()->MassForm(TP,__F,*GetFem(),__U);
    PatchDiscretization::LocalToGlobal(f,__F,iq,s);
  }
}

/**********************************************************/

void Q4::ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const
{
  int ncomp = u.ncomp();
  err.ncomp() = ncomp;
  err.reservesize(3);
  err = 0.;

  CompVector<double> lerr(ncomp,3);

  nmatrix<double> T;

  GlobalToGlobalData();
  ES->SetParameterData(__QP);

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); iq++)
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

/**********************************************************/

void Q4::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  RHS.SetParameterData(__QP);

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocalData(iq);
    GetIntegrator()->Rhs(RHS,__F,*GetFem(),__QN,__QC);
    PatchDiscretization::LocalToGlobal(f,__F,iq,s);
  }
}

/**********************************************************/

void Q4::BoundaryRhs(GlobalVector& f, const IntSet& Colors,  const BoundaryRightHandSide& NRHS, double s) const
{  
  int dim = GetPatchMesh()->dimension();
  nmatrix<double> T;

  GlobalToGlobalData();
  NRHS.SetParameterData(__QP);
  /// die cell2patch - liste

  const nvector<IntVector>& q4patch2cell = GetGascoigneMesh()->GetPatchIndexHandler().GetAllQ4Patch2Cell();

  nvector<int> cell2q4patch(GetMesh()->ncells());
  for(int p=0; p<q4patch2cell.size(); p++)
  {
    for(int i=0; i<q4patch2cell[p].size(); i++)
    {
      cell2q4patch[q4patch2cell[p][i]] = p;
    }
  }

  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;
    HASHSET<int> habschon;

    const IntVector& q = *GetMesh()->PatchOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
    for(int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2q4patch[iq];
      int ile = l[i];

      // gabs den patch schon?
      if(habschon.find((ip<<dim)+ile)!=habschon.end())
      {
        continue;
      }
      habschon.insert((ip<<dim)+ile);

      Transformation(T,ip);
      GetFem()->ReInit(T);

      GlobalToLocalData(ip);
      GetIntegrator()->BoundaryRhs(NRHS,__F,*GetFem(),ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,s);
    }
  }
}

/**********************************************************/

void Q4::InitFilter(nvector<double>& F) const
{
  PressureFilter* PF = static_cast<PressureFilter*>(&F);
  assert(PF);

  if (!PF->Active()) return;

  PF->ReInit(GetMesh()->nnodes(),GetMesh()->nhanging());
  nmatrix<double> T;

  int nv = GetPatchMesh()->nodes_per_q4patch();
  EntryMatrix  E(nv,1);

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    double cellsize = GetIntegrator()->MassMatrix(E,*GetFem());
    PF->AddDomainPiece(cellsize);

    nvector<int> ind = *GetPatchMesh()->IndicesOfQ4Patch(iq);
    HN->CondenseHanging(E,ind);

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

/**********************************************************/

double Q4::ComputeBoundaryFunctional(const GlobalVector& u, const IntSet& Colors, const BoundaryFunctional& BF) const
{
  int dim = GetPatchMesh()->dimension();
  nmatrix<double> T;

  GlobalToGlobalData();
  BF.SetParameterData(__QP);
  /// die cell2patch - liste

  const nvector<IntVector>& q4patch2cell = GetGascoigneMesh()->GetPatchIndexHandler().GetAllQ4Patch2Cell();

  nvector<int> cell2q4patch(GetMesh()->ncells());
  for(int p=0; p<q4patch2cell.size(); p++)
  {
    for(int i=0; i<q4patch2cell[p].size(); i++)
    {
      cell2q4patch[q4patch2cell[p][i]] = p;
    }
  }

  double j=0.;
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;
    HASHSET<int> habschon;

    const IntVector& q = *GetMesh()->PatchOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
    for(int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2q4patch[iq];
      int ile = l[i];

      // gabs den patch schon?
      if(habschon.find((ip<<dim)+ile)!=habschon.end())
      {
        continue;
      }
      habschon.insert((ip<<dim)+ile);

      Transformation(T,ip);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,ip);
      j += GetIntegrator()->ComputeBoundaryFunctional(BF,*GetFem(),ile,col,__U,__QN,__QC);
    }
  }
  return j;
}

/**********************************************************/

double Q4::ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const
{
  GlobalToGlobalData();
  F.SetParameterData(__QP);

  nmatrix<double> T;
  double j=0.;
  for(int iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocal(__U,u,iq);
    j += GetIntegrator()->ComputeDomainFunctional(F,*GetFem(),__U,__QN,__QC);
  }
  return j;
}

/**********************************************************/
}

#undef HASHSET

