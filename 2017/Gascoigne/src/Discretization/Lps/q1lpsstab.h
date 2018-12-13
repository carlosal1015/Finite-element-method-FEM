/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#ifndef  __Q1LpsStab_h
#define  __Q1LpsStab_h

#include  "patchdiscretization.h"
#include  "hnstructureinterface.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1LpsStab

////
////
/////////////////////////////////////////////

/*----------------------------------------------*/

class Q1LpsStab : public PatchDiscretization
{
 protected:

  const HNStructureInterface* HN;

  nvector<int> GetLocalIndices(int iq) const {
    return *GetPatchMesh()->IndicesOfPatch(iq);
  }
  void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

 public:

  Q1LpsStab() : PatchDiscretization() {};
  int n() const { return GetMesh()->nnodes();}
  int nc() const { return GetMesh()->ncells();}
  virtual void BasicInit(const ParamFile* paramfile, const HNStructureInterface*);
};

/*----------------------------------------------*/

class Q1LpsStab2d : public Q1LpsStab
{
 protected:

  void Transformation(FemInterface::Matrix& T, int iq) const;

 public:

  Q1LpsStab2d() : Q1LpsStab() {};
  void BasicInit(const ParamFile* paramfile, const HNStructureInterface*);
};

/*----------------------------------------------*/

class Q1LpsStab3d : public Q1LpsStab
{
 protected:

  void Transformation(FemInterface::Matrix& T, int iq) const;

 public:

  Q1LpsStab3d() : Q1LpsStab() {};
  void BasicInit(const ParamFile* paramfile, const HNStructureInterface*);
};

/*----------------------------------------------*/
}
#endif
