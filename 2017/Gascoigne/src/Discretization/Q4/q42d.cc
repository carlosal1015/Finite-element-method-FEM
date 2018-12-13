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


#include "q42d.h"
#include "baseq42d.h"
#include "finiteelement.h"
#include "galerkinintegratorq4.h"
#include "hnstructureq42d.h"
#include "transformation2d.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

Q42d::Q42d() : Q4()
{
  assert(HN==NULL);
  HN = new HNStructureQ42d;
}

/**********************************************************/

Q42d::~Q42d()
{
  if(HN)
  {
    delete HN;
    HN = NULL;
  }
}

/**********************************************************/

int Q42d::GetPatchNumber(const Vertex2d& p0, Vertex2d& p) const
{
  int iq;

  for(iq=0; iq<GetPatchMesh()->nq4patches(); ++iq)
  {
    bool found = true;
    const IntVector& IOP = GetPatchMesh()->CoarseIndicesQ4(iq);

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

void Q42d::VertexTransformation(const Vertex2d& p0, Vertex2d& p, int iq) const
{
  nmatrix<double> T;
  Transformation(T,iq);

  Transformation2d<BaseQ42d> Tr;
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

/**********************************************************/

void Q42d::BasicInit(const ParamFile* paramfile)
{
  if(GetIntegrator()==NULL)
  {
    PatchDiscretization::GetIntegratorPointer() = new GalerkinIntegratorQ4<2>;
  }
  assert(GetIntegrator());

  GetIntegratorPointer()->BasicInit();

  typedef Transformation2d<BaseQ42d>           TransQ4;
  typedef FiniteElement<2,1,TransQ4,BaseQ42d>  FiniteElement;

  if(GetFem()==NULL)
  {
    PatchDiscretization::GetFemPointer() = new FiniteElement;
  }
  assert(GetFem());

  PatchDiscretization::BasicInit(paramfile);
}

/**********************************************************/
}

