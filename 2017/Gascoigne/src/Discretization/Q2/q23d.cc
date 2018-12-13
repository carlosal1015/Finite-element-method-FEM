/**
*
* Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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


#include  "q23d.h"
#include  "galerkinintegratorq2.h"
#include  "transformation3d.h"
#include  "finiteelement.h"
#include  "baseq23d.h"
#include  "sparsestructure.h"
#include  "gascoignemeshtransfer.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "hnstructureq23d.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
Q23d::Q23d() : Q2()
{
  HN = new HNStructureQ23d;
}

/* ----------------------------------------- */

Q23d::~Q23d()
{
  if (HN) delete HN;
  HN = NULL;
}

/* ----------------------------------------- */

nmatrix<double> Q23d::GetLocalInterpolationWeights(int iq) const
{
  int nn = GetMesh()->nodes_per_cell(iq);
  nmatrix<double> w(nn,nn);
  w.zero();
  w(0,1) =  0.5 ; w(0,2) = 0.5 ; w(0,3) = 0.25; w(0,4) = 0.5  ; w(0,5) = 0.25 ; w(0,6) = 0.25 ; w(0,7) = 0.125;
  w(1,0) =  0.5 ; w(1,2) = 0.25; w(1,3) = 0.5 ; w(1,4) = 0.25 ; w(1,5) = 0.5  ; w(1,6) = 0.125; w(1,7) = 0.25;
  w(2,0) =  0.5 ; w(2,1) = 0.25; w(2,3) = 0.5 ; w(2,4) = 0.25 ; w(2,5) = 0.125; w(2,6) = 0.5  ; w(2,7) = 0.25;
  w(3,0) =  0.25; w(3,1) = 0.5 ; w(3,2) = 0.5 ; w(3,4) = 0.125; w(3,5) = 0.25 ; w(3,6) = 0.25 ; w(3,7) = 0.5;
  w(4,0) =  0.5 ; w(4,1) = 0.25; w(4,2) = 0.25; w(4,3) = 0.125; w(4,5) = 0.5  ; w(4,6) = 0.5  ; w(4,7) = 0.25;
  w(5,0) =  0.25; w(5,1) = 0.5 ; w(5,2) = 0.125;w(5,3) = 0.25 ; w(5,4) = 0.5  ; w(5,6) = 0.25 ; w(5,7) = 0.5;
  w(6,0) =  0.25; w(6,1) = 0.125;w(6,2) = 0.5 ; w(6,3) = 0.25 ; w(6,4) = 0.5  ; w(6,5) = 0.25 ; w(6,7) = 0.5;
  w(7,0) =  0.125;w(7,1) = 0.25; w(7,2) = 0.25; w(7,3) = 0.5  ; w(7,4) = 0.25 ; w(7,5) = 0.5  ; w(7,6) = 0.5;
  return w;
}

/* ----------------------------------------- */

int Q23d::GetPatchNumber(const Vertex3d& p0, Vertex3d& p) const
{
  int iq;
  
  for(iq=0; iq<GetPatchMesh()->npatches(); ++iq)
  {
    bool found = true;
    const IntVector& IOP = GetPatchMesh()->CoarseIndices(iq);
    
    for(int d=0; d<3; ++d)
    {
      double min=GetMesh()->vertex3d(IOP[0])[d];
      double max=min;
      for(int j=1; j<8; ++j)
      {
        double x = GetMesh()->vertex3d(IOP[j])[d];
        
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
    
    for(int d=0; d<3; ++d)
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

void Q23d::VertexTransformation(const Vertex3d& p0, Vertex3d& p, int iq) const
{
  nmatrix<double> T;
  Transformation(T,iq);

  Transformation3d<BaseQ23d> Tr;
  Tr.init(T);

  Vertex3d res;
  
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

void Q23d::BasicInit(const ParamFile* paramfile)
{
  if(!PatchDiscretization::GetIntegrator())
    PatchDiscretization::GetIntegratorPointer() =  new GalerkinIntegratorQ2<3>;
  assert(GetIntegrator());

  GetIntegratorPointer()->BasicInit();  

  if(!PatchDiscretization::GetFem())
    {
      typedef Transformation3d<BaseQ23d>           TransQ2;
      typedef FiniteElement<3,2,TransQ2,BaseQ23d>  FiniteElement;
      
      PatchDiscretization::GetFemPointer() =  new FiniteElement;
    }
  assert(GetFem());
  
  PatchDiscretization::BasicInit(paramfile);
}

/* ----------------------------------------- */

void Q23d::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
{
  MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);

  assert(IP);
  IP->BasicInit(MT);
  return;
}
}
