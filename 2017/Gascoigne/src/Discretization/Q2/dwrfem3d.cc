/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2008 by the Gascoigne 3D authors
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


#include "dwrfem.h" 
#include "galerkinintegratorq2.h"
#include "baseq23d.h"
#include "integratorq1q2.h"
#include "hnstructureq23d.h"

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

/*---------------------------------------------------*/

DwrFem3d::DwrFem3d() : Q23d()
{
  HNLow = new HNStructureQ13d;
}

/*---------------------------------------------------*/

DwrFem3d::~DwrFem3d()
{
  if(HNLow)
  {
    delete HNLow;
    HNLow = NULL;
  }
}

/*---------------------------------------------------*/

void DwrFem3d::BasicInit(const ParamFile* paramfile)
{
  assert(PatchDiscretization::GetIntegrator()==NULL);
  PatchDiscretization::GetIntegratorPointer() = new IntegratorQ1Q2<3>;

  GetIntegratorPointer()->BasicInit();

  assert(PatchDiscretization::GetFem()==NULL);
  typedef Transformation3d<BaseQ23d>          TransQ2;
  typedef FiniteElement<3,2,TransQ2,BaseQ23d> FiniteElement;

  PatchDiscretization::GetFemPointer() = new FiniteElement;
  PatchDiscretization::BasicInit(paramfile);
}

/*---------------------------------------------------*/

void DwrFem3d::ReInit(const MeshInterface* MP)
{
  Q23d::ReInit(MP);

  HNLow->ReInit(MP);
}

/*---------------------------------------------------*/

void DwrFem3d::TransformationQ1(FemInterface::Matrix& T, int iq) const
{
  int dim = GetMesh()->dimension();
  int ne = GetMesh()->nodes_per_cell(iq);

  nvector<int> indices = GetPatchMesh()->CoarseIndices(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);
  for(int ii=0;ii<ne;ii++)
    {
      Vertex3d v = GetMesh()->vertex3d(indices[ii]);
      T(0,ii) = v.x();               
      T(1,ii) = v.y();
      T(2,ii) = v.z();
    }
}

/*---------------------------------------------------*/
/*---------------------------------------------------*/

void DwrFemQ1Q23d::DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex3d& p0, int i, double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  Vertex3d Tranfo_p0;
   
  int iq = GetPatchNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "DwrFemQ1Q23d::DiracRhsPoint point not found\n";
      abort();
    }

  nmatrix<double> TH,TL;

  Transformation  (TH,iq);
  TransformationQ1(TL,iq);

  const FemInterface& HighOrderFem(*GetFem());

  HighOrderFem.ReInit(TH);
  LowOrderFem .ReInit(TL);
  
  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);
  
  GlobalToLocalData(iq);
  GlobalToGlobalData();
  DRHS.SetParameterData(__QP);

  I->DiracRhsPoint(__F,HighOrderFem,LowOrderFem,Tranfo_p0,DRHS,i,__QN,__QC);
  PatchDiscretization::LocalToGlobal(f,__F,iq,s);
}

/*---------------------------------------------------*/

void DwrFemQ1Q23d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,iq);
      I->Form(EQ,__F,HighOrderFem,LowOrderFem,__U,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void DwrFemQ1Q23d::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,iq);
      I->AdjointForm(EQ,__F,HighOrderFem,LowOrderFem,__U,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void DwrFemQ1Q23d::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, 
    const BoundaryEquation& BE, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  BE.SetParameterData(__QP);
  
  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;
    
    const IntVector& q = *GetMesh()->PatchOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
    for (int i=0; i<q.size(); i++)
    {
      int ip  = q[i];
      int ile = l[i];

      Transformation  (TH,ip);
      TransformationQ1(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,ip);
      I->BoundaryForm(BE,__F,HighOrderFem,LowOrderFem,__U,ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,d);
    }
  }
}

/* ----------------------------------------- */

void DwrFemQ1Q23d::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  RHS.SetParameterData(__QP);

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocalData(iq);
      I->Rhs(RHS,__F,HighOrderFem,LowOrderFem,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,iq,s);
    }
}

/* ----------------------------------------- */

void DwrFemQ1Q23d::BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  NRHS.SetParameterData(__QP);
  
  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;

    const IntVector& q = *GetMesh()->PatchOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
    for (int i=0; i<q.size(); i++)
    {
      int ip = q[i];
      int ile = l[i];

      Transformation  (TH,ip);
      TransformationQ1(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocalData(ip);
      I->BoundaryRhs(NRHS,__F,HighOrderFem,LowOrderFem,ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,s);
    }
  }
}

/* ----------------------------------------- */

void DwrFemQ1Q23d::MassMatrix(MatrixInterface& A) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      I->MassMatrix(__E,HighOrderFem,LowOrderFem);
      LocalToGlobal(A,__E,iq,1.);
    }
}

/* ----------------------------------------- */

void Gascoigne::DwrFemQ1Q23d::MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,iq);
      I->MassForm(TP,__F,HighOrderFem,LowOrderFem,__U);
      PatchDiscretization::LocalToGlobal(f,__F,iq,s);
    }
}

/* ----------------------------------------- */

void DwrFemQ1Q23d::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  dynamic_cast<HNStructureQ23d*>(HN)->CondenseHangingLowerHigher(E,indices);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,E,s);
}

/*---------------------------------------------------*/
/*---------------------------------------------------*/

void DwrFemQ2Q13d::DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex3d& p0, int i, double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  Vertex3d Tranfo_p0;
   
  int iq = GetPatchNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "DwrFemQ2Q13d::DiracRhsPoint point not found\n";
      abort();
    }

  nmatrix<double> TH,TL;

  Transformation  (TH,iq);
  TransformationQ1(TL,iq);

  const FemInterface& HighOrderFem(*GetFem());

  HighOrderFem.ReInit(TH);
  LowOrderFem .ReInit(TL);
  
  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);
  
  GlobalToLocalData(iq);
  GlobalToGlobalData();
  DRHS.SetParameterData(__QP);

  I->DiracRhsPoint(__F,LowOrderFem,HighOrderFem,Tranfo_p0,DRHS,i,__QN,__QC);
  PatchDiscretization::LocalToGlobal(f,__F,iq,s);
}

/*---------------------------------------------------*/

void DwrFemQ2Q13d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,iq);
      I->Form(EQ,__F,LowOrderFem,HighOrderFem,__U,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void DwrFemQ2Q13d::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__QP);

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,iq);
      I->AdjointForm(EQ,__F,LowOrderFem,HighOrderFem,__U,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void DwrFemQ2Q13d::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, 
    const BoundaryEquation& BE, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  BE.SetParameterData(__QP);
  
  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  const nvector<nvector<int> >& patch2cell  =
      GetGascoigneMesh()->GetPatchIndexHandler().GetAllPatch2Cell();

  nvector<int> cell2patch(GetMesh()->ncells());
  for (int p=0;p<patch2cell.size();++p)
    for (int i=0;i<patch2cell[p].size();++i)
       cell2patch[patch2cell[p][i]]=p;
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;

    HASHSET<int> habschon;
    
    const IntVector& q = *GetMesh()->CellOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalOnBoundary(col);
    for (int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2patch[iq];
      
      // gabs den patch schon?
      if (habschon.find(ip)!=habschon.end()) continue;
      habschon.insert(ip);

      int ile = l[i];

      Transformation  (TH,ip);
      TransformationQ1(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,ip);
      I->BoundaryForm(BE,__F,LowOrderFem,HighOrderFem,__U,ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,d);
    }
  }
}

/* ----------------------------------------- */

void DwrFemQ2Q13d::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  RHS.SetParameterData(__QP);

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocalData(iq);
      I->Rhs(RHS,__F,LowOrderFem,HighOrderFem,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,iq,s);
    }
}

/* ----------------------------------------- */

void DwrFemQ2Q13d::BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  NRHS.SetParameterData(__QP);
  
  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  const nvector<nvector<int> >& patch2cell  =
      GetGascoigneMesh()->GetPatchIndexHandler().GetAllPatch2Cell();

  nvector<int> cell2patch(GetMesh()->ncells());
  for (int p=0;p<patch2cell.size();++p)
    for (int i=0;i<patch2cell[p].size();++i)
       cell2patch[patch2cell[p][i]]=p;
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;
    HASHSET<int> habschon;

    const IntVector& q = *GetMesh()->CellOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalOnBoundary(col);
    for (int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2patch[iq];
      
      // gabs den patch schon?
      if (habschon.find(ip)!=habschon.end()) continue;
      habschon.insert(ip);
      
      int ile = l[i];

      Transformation  (TH,ip);
      TransformationQ1(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocalData(ip);
      I->BoundaryRhs(NRHS,__F,LowOrderFem,HighOrderFem,ile,col,__QN,__QC);
      PatchDiscretization::LocalToGlobal(f,__F,ip,s);
    }
  }
}

/* ----------------------------------------- */

void DwrFemQ2Q13d::MassMatrix(MatrixInterface& A) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      I->MassMatrix(__E,LowOrderFem,HighOrderFem);
      LocalToGlobal(A,__E,iq,1.);
    }
}

/* ----------------------------------------- */

void Gascoigne::DwrFemQ2Q13d::MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,iq);
      I->MassForm(TP,__F,LowOrderFem,HighOrderFem,__U);
      PatchDiscretization::LocalToGlobal(f,__F,iq,s);
    }
}

/* ----------------------------------------- */

void DwrFemQ2Q13d::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  dynamic_cast<HNStructureQ23d*>(HN)->CondenseHangingHigherLower(E,indices);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,E,s);
}

}
#undef HASHSET
