/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#include  "q1lpsstab.h"
#include  "transformation2d.h"
#include  "transformation3d.h"
#include  "finiteelement.h"
#include  "baseq1patch.h"
#include  "baseq13dpatch.h"
#include  "lpsintegrator.h"

namespace Gascoigne
{
/* ----------------------------------------- */

void Q1LpsStab::BasicInit(const ParamFile* paramfile, const HNStructureInterface* hn)
{
  PatchDiscretization::BasicInit(paramfile);
  HN = hn;
  assert(HN);
}

/* ----------------------------------------- */

void Q1LpsStab::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  HN->CondenseHangingPatch(E,indices);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,__E,s);
}

/* ----------------------------------------- */

void Q1LpsStab2d::BasicInit(const ParamFile* paramfile, const HNStructureInterface* hn)
{
  assert(PatchDiscretization::GetIntegrator()==NULL);
  PatchDiscretization::GetIntegratorPointer() =  new LpsIntegratorQ1<2>;
  
  assert(PatchDiscretization::GetFem()==NULL);
  typedef Transformation2d<BaseQ12d>           TransQ1;
  typedef FiniteElement<2,1,TransQ1,BaseQ12dPatch>  FiniteElement;
  PatchDiscretization::GetFemPointer() =  new FiniteElement;
  
  Q1LpsStab::BasicInit(paramfile,hn);  
}

/* ----------------------------------------- */

void Q1LpsStab3d::BasicInit(const ParamFile* paramfile, const HNStructureInterface* hn)
{
  assert(PatchDiscretization::GetIntegrator()==NULL);
  PatchDiscretization::GetIntegratorPointer() =  new LpsIntegratorQ1<3>;
  
  assert(PatchDiscretization::GetFem()==NULL);
  typedef Transformation3d<BaseQ13d>           TransQ1;
  typedef FiniteElement<3,2,TransQ1,BaseQ13dPatch>  FiniteElement;
  PatchDiscretization::GetFemPointer() =  new FiniteElement;
  
  Q1LpsStab::BasicInit(paramfile,hn);  
}

/* ----------------------------------------- */

void Q1LpsStab2d::Transformation(FemInterface::Matrix& T, int iq) const
{
  int dim = 2;
  int ne = GetPatchMesh()->nodes_per_cell(iq);

  nvector<int> indices = GetPatchMesh()->CoarseIndices(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);

  for(int ii=0;ii<ne;ii++)
    {
      Vertex2d v = GetPatchMesh()->vertex2d(indices[ii]);
      T(0,ii) = v.x();               
      T(1,ii) = v.y();
    }
}

/* ----------------------------------------- */

void Q1LpsStab3d::Transformation(FemInterface::Matrix& T, int iq) const
{
  int dim = 3;
  int ne = GetPatchMesh()->nodes_per_cell(iq);

  nvector<int> indices = GetPatchMesh()->CoarseIndices(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);

  for(int ii=0;ii<ne;ii++)
    {
      Vertex3d v = GetPatchMesh()->vertex3d(indices[ii]);
      T(0,ii) = v.x();               
      T(1,ii) = v.y();
      T(2,ii) = v.z();
    }
}

}
