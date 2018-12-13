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


#include "dwrfemq2.h"
#include "baseq42d.h"
#include "hnstructureq42d.h"
#include "integratorq2q4.h"
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

DwrFemQ22d::DwrFemQ22d() : Q42d()
{
  HNLow = new HNStructureQ22d;
}

/**********************************************************/

DwrFemQ22d::~DwrFemQ22d()
{
  if(HNLow)
  {
    delete HNLow;
    HNLow = NULL;
  }
}

/**********************************************************/

void DwrFemQ22d::TransformationQ2(FemInterface::Matrix& T, int iq) const
{
  int dim = GetMesh()->dimension();
  assert(dim==2);
  int ne = GetPatchMesh()->nodes_per_patch();

  IntVector indices = GetPatchMesh()->Q2IndicesOfQ4Patch(iq);

  assert(ne==indices.size());

  T.memory(dim,ne);
  for(int ii=0; ii<ne; ii++)
  {
    Vertex2d v = GetMesh()->vertex2d(indices[ii]);
    T(0,ii) = v.x();
    T(1,ii) = v.y();
  }
}

/**********************************************************/

void DwrFemQ22d::BasicInit(const ParamFile* paramfile)
{
  assert(PatchDiscretization::GetIntegrator()==NULL);
  PatchDiscretization::GetIntegratorPointer() = new IntegratorQ2Q4<2>;
  assert(PatchDiscretization::GetIntegrator());

  GetIntegratorPointer()->BasicInit();

  assert(PatchDiscretization::GetFem()==NULL);
  typedef Transformation2d<BaseQ42d>          TransQ4;
  typedef FiniteElement<2,1,TransQ4,BaseQ42d> FiniteElement;

  PatchDiscretization::GetFemPointer() = new FiniteElement;
  PatchDiscretization::BasicInit(paramfile);
}

/**********************************************************/

void DwrFemQ22d::ReInit(const MeshInterface* MP)
{
  Q42d::ReInit(MP);

  HNLow->ReInit(MP);
}

/**********************************************************/
/**********************************************************/

void DwrFemQ2Q42d::DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex2d& p0, int i, double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  Vertex2d Tranfo_p0;

  int iq = GetPatchNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "DwrFemQ2Q42d::DiracRhsPoint point not found\n";
      abort();
    }

  nmatrix<double> TH,TL;

  Transformation  (TH,iq);
  TransformationQ2(TL,iq);

  const FemInterface& HighOrderFem(*GetFem());

  HighOrderFem.ReInit(TH);
  LowOrderFem .ReInit(TL);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  GlobalToLocalData(iq);
  GlobalToGlobalData();
  DRHS.SetParameterData(__QP);

  I->DiracRhsPoint(__F,HighOrderFem,LowOrderFem,Tranfo_p0,DRHS,i,__QN,__QC);
  PatchDiscretization::LocalToGlobal(f,__F,iq,s);
}

/**********************************************************/

void DwrFemQ2Q42d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->nq4patches();++iq)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    GlobalToLocal(__U,u,iq);
    I->Form(EQ,__F,HighOrderFem,LowOrderFem,__U,__QN,__QC);
    PatchDiscretization::LocalToGlobal(f,__F,iq,d);
  }
}

/**********************************************************/

void DwrFemQ2Q42d::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->nq4patches();++iq)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    GlobalToLocal(__U,u,iq);
    I->AdjointForm(EQ,__F,HighOrderFem,LowOrderFem,__U,__QN,__QC);
    PatchDiscretization::LocalToGlobal(f,__F,iq,d);
  }
}

/**********************************************************/

void DwrFemQ2Q42d::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  BE.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  const nvector<IntVector>& q4patch2cell = GetGascoigneMesh()->GetPatchIndexHandler().GetAllQ4Patch2Cell();

  nvector<int> cell2q4patch(GetMesh()->ncells());
  for(int p=0; p<q4patch2cell.size(); p++)
  {
    for(int i=0; i<q4patch2cell[p].size(); i++)
    {
      cell2q4patch[q4patch2cell[p][i]] = p;
    }
  }

  for(IntSet::const_iterator p=Colors.begin(); p!=Colors.end(); p++)
  {
    int col = *p;

    HASHSET<int> habschon;

    const IntVector& q = *GetMesh()->CellOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalOnBoundary(col);
    for(int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2q4patch[iq];
      int ile = l[i];

      // gabs den patch schon?
      if(habschon.find((ip<<2)+ile)!=habschon.end())
      {
        continue;
      }
      habschon.insert((ip<<2)+ile);

      Transformation  (TH,ip);
      TransformationQ2(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,ip);
      I->BoundaryForm(BE,__F,HighOrderFem,LowOrderFem,__U,ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,d);
    }
  }
}

/**********************************************************/

void DwrFemQ2Q42d::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  RHS.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); iq++)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    GlobalToLocalData(iq);
    I->Rhs(RHS,__F,HighOrderFem,LowOrderFem,__QN,__QC);
    PatchDiscretization::LocalToGlobal(f,__F,iq,s);
  }
}

/**********************************************************/

void DwrFemQ2Q42d::BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  NRHS.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  const nvector<IntVector>& q4patch2cell = GetGascoigneMesh()->GetPatchIndexHandler().GetAllQ4Patch2Cell();

  IntVector cell2q4patch(GetMesh()->ncells());
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

    const IntVector& q = *GetMesh()->CellOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalOnBoundary(col);
    for(int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2q4patch[iq];
      int ile = l[i];

      // gabs den patch schon?
      if(habschon.find((ip<<2)+ile)!=habschon.end())
      {
        continue;
      }
      habschon.insert((ip<<2)+ile);

      Transformation  (TH,ip);
      TransformationQ2(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocalData(ip);
      I->BoundaryRhs(NRHS,__F,HighOrderFem,LowOrderFem,ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,s);
    }
  }
}

/**********************************************************/

void DwrFemQ2Q42d::MassMatrix(MatrixInterface& M) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    I->MassMatrix(__E,HighOrderFem,LowOrderFem);
    LocalToGlobal(M,__E,iq,1.);
  }
}

/**********************************************************/

void DwrFemQ2Q42d::MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->nq4patches();++iq)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    GlobalToLocal(__U,u,iq);
    I->MassForm(TP,__F,HighOrderFem,LowOrderFem,__U);
    PatchDiscretization::LocalToGlobal(f,__F,iq,s);
  }
}

/**********************************************************/

void DwrFemQ2Q42d::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  dynamic_cast<HNStructureQ42d*>(HN)->CondenseHangingLowerHigher(E,indices);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,E,s);
}

/**********************************************************/
/**********************************************************/

void DwrFemQ4Q22d::DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex2d& p0, int i, double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  Vertex2d Tranfo_p0;

  int iq = GetPatchNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "DwrFemQ4Q22d::DiracRhsPoint point not found\n";
      abort();
    }

  nmatrix<double> TH,TL;

  Transformation  (TH,iq);
  TransformationQ2(TL,iq);

  const FemInterface& HighOrderFem(*GetFem());

  HighOrderFem.ReInit(TH);
  LowOrderFem .ReInit(TL);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  GlobalToLocalData(iq);
  GlobalToGlobalData();
  DRHS.SetParameterData(__QP);

  I->DiracRhsPoint(__F,LowOrderFem,HighOrderFem,Tranfo_p0,DRHS,i,__QN,__QC);
  PatchDiscretization::LocalToGlobal(f,__F,iq,s);
}

/**********************************************************/

void DwrFemQ4Q22d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->nq4patches();++iq)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    GlobalToLocal(__U,u,iq);
    I->Form(EQ,__F,LowOrderFem,HighOrderFem,__U,__QN,__QC);
    PatchDiscretization::LocalToGlobal(f,__F,iq,d);
  }
}

/**********************************************************/

void DwrFemQ4Q22d::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->nq4patches();++iq)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    GlobalToLocal(__U,u,iq);
    I->AdjointForm(EQ,__F,LowOrderFem,HighOrderFem,__U,__QN,__QC);
    PatchDiscretization::LocalToGlobal(f,__F,iq,d);
  }
}

/**********************************************************/

void DwrFemQ4Q22d::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  BE.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  const nvector<IntVector>& q4patch2cell = GetGascoigneMesh()->GetPatchIndexHandler().GetAllQ4Patch2Cell();

  nvector<int> cell2q4patch(GetMesh()->ncells());
  for(int p=0; p<q4patch2cell.size(); p++)
  {
    for(int i=0; i<q4patch2cell[p].size(); i++)
    {
      cell2q4patch[q4patch2cell[p][i]] = p;
    }
  }

  for(IntSet::const_iterator p=Colors.begin(); p!=Colors.end(); p++)
  {
    int col = *p;

    HASHSET<int> habschon;

    const IntVector& q = *GetMesh()->CellOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalOnBoundary(col);
    for(int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2q4patch[iq];

      // gabs den patch schon?
      if(habschon.find(ip)!=habschon.end())
      {
        continue;
      }
      habschon.insert(ip);

      int ile = l[i];

      Transformation  (TH,ip);
      TransformationQ2(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,ip);
      I->BoundaryForm(BE,__F,LowOrderFem,HighOrderFem,__U,ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,d);
    }
  }
}

/**********************************************************/

void DwrFemQ4Q22d::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  RHS.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); iq++)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    GlobalToLocalData(iq);
    I->Rhs(RHS,__F,LowOrderFem,HighOrderFem,__QN,__QC);
    PatchDiscretization::LocalToGlobal(f,__F,iq,s);
  }
}

/**********************************************************/

void DwrFemQ4Q22d::BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  NRHS.SetParameterData(__QP);

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  const nvector<IntVector>& q4patch2cell = GetGascoigneMesh()->GetPatchIndexHandler().GetAllQ4Patch2Cell();

  IntVector cell2q4patch(GetMesh()->ncells());
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

    const IntVector& q = *GetMesh()->CellOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalOnBoundary(col);
    for(int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2q4patch[iq];

      // gabs den patch schon?
      if(habschon.find(ip)!=habschon.end())
      {
        continue;
      }
      habschon.insert(ip);

      int ile = l[i];

      Transformation  (TH,ip);
      TransformationQ2(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocalData(ip);
      I->BoundaryRhs(NRHS,__F,LowOrderFem,HighOrderFem,ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,s);
    }
  }
}

/**********************************************************/

void DwrFemQ4Q22d::MassMatrix(MatrixInterface& M) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    I->MassMatrix(__E,LowOrderFem,HighOrderFem);
    LocalToGlobal(M,__E,iq,1.);
  }
}

/**********************************************************/

void DwrFemQ4Q22d::MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ2Q4<2>* I = dynamic_cast<const IntegratorQ2Q4<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->nq4patches();++iq)
  {
    Transformation  (TH,iq);
    TransformationQ2(TL,iq);

    HighOrderFem.ReInit(TH);
    LowOrderFem .ReInit(TL);

    GlobalToLocal(__U,u,iq);
    I->MassForm(TP,__F,LowOrderFem,HighOrderFem,__U);
    PatchDiscretization::LocalToGlobal(f,__F,iq,s);
  }
}

/**********************************************************/

void DwrFemQ4Q22d::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  dynamic_cast<HNStructureQ42d*>(HN)->CondenseHangingHigherLower(E,indices);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,E,s);
}

/**********************************************************/
}

#undef HASHSET
