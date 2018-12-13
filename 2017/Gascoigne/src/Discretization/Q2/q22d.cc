/**
*
* Copyright (C) 2004, 2006 by the Gascoigne 3D authors
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


#include  "q22d.h"
#include  "galerkinintegratorq2.h"
#include  "transformation2d.h"
#include  "finiteelement.h"
#include  "baseq22d.h"
#include  "sparsestructure.h"
#include  "gascoignemeshtransfer.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "hnstructureq22d.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
Q22d::Q22d() : Q2()
{
  assert(HN==NULL);
  HN = new HNStructureQ22d;
}

/* ----------------------------------------- */

Q22d::~Q22d()
{
  if (HN) delete HN;
  HN = NULL;
}

/* ----------------------------------------- */

nmatrix<double> Q22d::GetLocalInterpolationWeights(int iq) const
{
  int nn = GetMesh()->nodes_per_cell(iq);
  nmatrix<double> w(nn,nn);
  w.zero();
  w(0,1) =  0.5  ; w(0,2) =  0.5  ; w(0,3) =  0.25;
  w(1,0) =  0.5  ; w(1,2) =  0.25 ; w(1,3) =  0.5 ;
  w(2,0) =  0.5  ; w(2,1) =  0.25 ; w(2,3) =  0.5 ;
  w(3,0) =  0.25 ; w(3,1) =  0.5  ; w(3,2) =  0.5 ;
  return w;
}

/* ----------------------------------------- */

int Q22d::GetPatchNumber(const Vertex2d& p0, Vertex2d& p) const
{
  int iq;
  
  for(iq=0; iq<GetPatchMesh()->npatches(); ++iq)
  {
    bool found = true;
    const IntVector& IOP = GetPatchMesh()->CoarseIndices(iq);
    
    for(int d=0; d<2; ++d)
    {
      double min=GetMesh()->vertex2d(IOP[0])[d];
      double max=min;
      for(int j=1; j<4; ++j)
      {
        double x = GetMesh()->vertex2d(IOP[j])[d];
        
        min = Gascoigne::min(min,x);
        max = Gascoigne::max(max,x);
      }
      if((p0[d]<min)||(p0[d]>max)) 
      {
        found = false;
        break;
      }
    }
    
    if(!found)
    {
      continue;
    }
    
    VertexTransformation(p0,p,iq);
    
    for(int d=0; d<2; ++d)
    {
      if((p[d]<0.-1.e-12)||(p[d]>1.+1.e-12))
      {
        found = false;
      }
    }
    
    if(found)
    {
      break;
    }
  }

  if(iq<GetPatchMesh()->npatches())
  {
    return iq;
  }
  else
  {
    return -1;
  }
}

/* ----------------------------------------- */

void Q22d::VertexTransformation(const Vertex2d& p0, Vertex2d& p, int iq) const
{
  nmatrix<double> T;
  Transformation(T,iq);

  Transformation2d<BaseQ22d> Tr;
  Tr.init(T);

  Vertex2d res;
  
  p = 0.5;
  
  for(int niter=1; ;niter++)
  {
    Tr.point(p);
    
    res = p0;
    res.add(-1,Tr.x());

    if(res.norm()<1.e-13)
    {
      break;
    }
    assert(niter<10);
    
    Tr.DTI().mult_ad(p,res);
  } 
}

/* ----------------------------------------- */

void Q22d::BasicInit(const ParamFile* paramfile)
{
  if (GetIntegrator()==NULL)
    PatchDiscretization::GetIntegratorPointer() =  new GalerkinIntegratorQ2<2>;
  assert(GetIntegrator());

  GetIntegratorPointer()->BasicInit();

  typedef Transformation2d<BaseQ22d>           TransQ2;
  typedef FiniteElement<2,1,TransQ2,BaseQ22d>  FiniteElement;

  if (GetFem()==NULL)
    PatchDiscretization::GetFemPointer() =  new FiniteElement;
  assert(GetFem());

  PatchDiscretization::BasicInit(paramfile);
}

/* ----------------------------------------- */

void Q22d::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
{
  {
    MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);
    if(IP)
      {
	IP->BasicInit(MT);
	return;
      }
  }
  MgInterpolatorMatrix* IP = dynamic_cast<MgInterpolatorMatrix*>(I);
  assert(IP);
  const GascoigneMeshTransfer* GMT = dynamic_cast<const GascoigneMeshTransfer*>(MT);
  assert(GMT);
  //
  // dast ist einfach von Q12d kopiert !!!!!
  //

  const map<int,fixarray<2,int> >& zweier = GMT->GetZweier();
  const map<int,fixarray<4,int> >& vierer = GMT->GetVierer();
  const map<int,fixarray<8,int> >& achter = GMT->GetAchter();
  const nvector<int>& c2f    = GMT->GetC2f();

  int n  = c2f.size() +   zweier.size() +   vierer.size() +   achter.size();
  int nt = c2f.size() + 2*zweier.size() + 4*vierer.size() + 8*achter.size();

  ColumnStencil& ST = IP->GetStencil();
  nvector<double>& val = IP->GetAlpha();

  SparseStructure SS;

  SS.build_begin(n);
  for(int i=0;i<c2f.size();i++)
    {
      SS.build_add(i,c2f[i]);
    }
  for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      for(int ii=0;ii<2;ii++) SS.build_add(il,n2[ii]);
    }
  for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    for(int ii=0;ii<4;ii++) SS.build_add(il,n4[ii]);
  }
  for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) {
    int il = p->first;
    fixarray<8,int> n8 = p->second;
    for(int ii=0;ii<8;ii++) SS.build_add(il,n8[ii]);
  }
  SS.build_end();

  assert(nt==SS.ntotal());

  ST.memory(&SS);

  val.reservesize(nt);

  for(int i=0;i<c2f.size();i++)
    {
      val[ST.Find(c2f[i],i)] = 1.;
    }
  for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      val[ST.Find(il,n2[0])] = 0.5;
      val[ST.Find(il,n2[1])] = 0.5;
    }
  for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    val[ST.Find(il,n4[0])] = 0.25;
    val[ST.Find(il,n4[1])] = 0.25;
    val[ST.Find(il,n4[2])] = 0.25;
    val[ST.Find(il,n4[3])] = 0.25;
  }
  for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) {
    int il = p->first;
    fixarray<8,int> n8 = p->second;
    for (int i=0; i<8; i++)
      {
	val[ST.Find(il,n8[i])] = 0.125;
      }
  }
}
}
