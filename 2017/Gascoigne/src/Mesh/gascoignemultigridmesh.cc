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


#include  "gascoignemultigridmesh.h"
#include  "gascoignemesh2d.h"
#include  "gascoignemesh3d.h"
#include  "gascoignemeshtransfer2d.h"
#include  "gascoignemeshtransfer3d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
GascoigneMultiGridMesh::GascoigneMultiGridMesh()
{
}

/*-----------------------------------------*/

GascoigneMultiGridMesh::~GascoigneMultiGridMesh()
{
  for(int l=0;l<M.size();l++)  
    {
      if (M[l]!=NULL)  {delete M[l]; M[l]=NULL;}
    }
  for(int l=0;l<T.size();l++)  
    {
      if (T[l]!=NULL)  {delete T[l]; T[l]=NULL;}
    }
}

/*-----------------------------------------*/

GascoigneMesh* GascoigneMultiGridMesh::NewMesh(int dim)
{
  assert(2<=dim && dim<=3);

  if(dim==2)
    return new GascoigneMesh2d;
  else if(dim==3)
    return new GascoigneMesh3d;

  abort();
}

/*-----------------------------------------*/

GascoigneMeshTransfer* GascoigneMultiGridMesh::NewTransfer(int dim)
{
  assert(2<=dim && dim<=3);

  if(dim==2)
    {
      return new GascoigneMeshTransfer2d;
    }
  else if(dim==3)
    {
      return new GascoigneMeshTransfer3d;
    }

  abort();
}

/*-----------------------------------------*/

void GascoigneMultiGridMesh::ReInit(int dim, int nlevels)
{
  // Mesh
  for(int l=0;l<M.size();l++)
    {
      if (M[l]!=NULL) { delete M[l]; M[l]=NULL;}
    }
  M.resize(nlevels);
  for(int l=0;l<M.size();l++)
    {
      M[l] = NewMesh(dim);
    }
  // Transfer
  for(int l=0;l<T.size();l++)
    {
      if(T[l]!=NULL) {delete T[l]; T[l]=NULL;}
    }
  T.resize(nlevels-1);
  for(int l=0;l<T.size();l++)
    {
      T[l] = NewTransfer(dim);
    }
}
}
