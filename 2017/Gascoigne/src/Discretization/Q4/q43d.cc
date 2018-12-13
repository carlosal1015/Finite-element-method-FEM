/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#include "q43d.h"
#include "baseq43d.h"
#include "finiteelement.h"
#include "galerkinintegratorq4.h"
#include "hnstructureq43d.h"
#include "transformation3d.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

Q43d::Q43d() : Q4()
{
  assert(HN==NULL);
  HN = new HNStructureQ43d;
}

/**********************************************************/

Q43d::~Q43d()
{
  if(HN)
  {
    delete HN;
    HN = NULL;
  }
}

/**********************************************************/

int Q43d::GetPatchNumber(const Vertex3d& p0, Vertex3d& p) const
{
  int iq;

  for(iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    bool found = true;
    const IntVector& IOP = GetPatchMesh()->CoarseIndicesQ4(iq);

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

  if(iq<GetPatchMesh()->nq4patches())
  {
    return iq;
  }
  else
  {
    return -1;
  }
}

/**********************************************************/

void Q43d::VertexTransformation(const Vertex3d& p0, Vertex3d& p, int iq) const
{
  nmatrix<double> T;
  Transformation(T,iq);

  Transformation3d<BaseQ43d> Tr;
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

/**********************************************************/

void Q43d::BasicInit(const ParamFile* paramfile)
{
  if(GetIntegrator()==NULL)
  {
    PatchDiscretization::GetIntegratorPointer() = new GalerkinIntegratorQ4<3>;
  }
  assert(GetIntegrator());

  GetIntegratorPointer()->BasicInit();

  typedef Transformation3d<BaseQ43d>           TransQ4;
  typedef FiniteElement<3,2,TransQ4,BaseQ43d>  FiniteElement;

  if(GetFem()==NULL)
  {
    PatchDiscretization::GetFemPointer() = new FiniteElement;
  }
  assert(GetFem());

  PatchDiscretization::BasicInit(paramfile);
}

/**********************************************************/
}

