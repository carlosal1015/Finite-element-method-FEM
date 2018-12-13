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


#ifndef  __GascoigneMultiGridMesh_h
#define  __GascoigneMultiGridMesh_h

#include  <vector>
#include  "gascoignemesh.h"
#include  "gascoignemeshtransfer.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class GascoigneMultiGridMesh
{
protected:

  std::vector<GascoigneMesh*>          M;
  std::vector<GascoigneMeshTransfer*>  T;

  virtual GascoigneMesh* NewMesh(int dim);
  virtual GascoigneMeshTransfer* NewTransfer(int dim);

public:
  
  virtual void ReInit(int dim, int nlevels);

  GascoigneMultiGridMesh();
  virtual ~GascoigneMultiGridMesh();

  int nlevels() const {return M.size();}

  const MeshInterface& operator()(int l) const {
    assert((l>=0)&&(l<M.size()));
    return *M[l];}

  const GascoigneMesh* GetGascoigneMesh(int l) const {
    assert((l>=0)&&(l<M.size()));
    return M[l];}
  GascoigneMesh*       GetGascoigneMesh(int l) {
    assert((l>=0)&&(l<M.size()));
    return M[l];}

  GascoigneMeshTransfer* GetTransfer(int l) {
    assert((l>=0)&&(l<T.size()));
    return T[l];}
  const GascoigneMeshTransfer* GetTransfer(int l) const {
    assert((l>=0)&&(l<T.size()));
    return T[l];}
};
}

#endif
