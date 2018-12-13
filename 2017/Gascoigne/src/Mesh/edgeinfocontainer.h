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


#ifndef __EdgeInfoContainer_h
#define __EdgeInfoContainer_h

#include "edgeinfo.h"
#include "hierarchicalmesh2d.h"
#include "hierarchicalmesh3d.h"
#include "nvector.h"

/**********************************************************/

namespace Gascoigne
{

class EdgeInfoContainerInterface
{
 public:

  EdgeInfoContainerInterface() {}
  virtual ~EdgeInfoContainerInterface() {}

  virtual const HierarchicalMesh* GetMesh() const=0;
  virtual void BasicInit(const HierarchicalMesh*, int)=0;
};

template<int DIM>
  class EdgeInfoContainer : public virtual EdgeInfoContainerInterface, public nvector<EdgeInfo<DIM>*>
{

 protected:

  const HierarchicalMesh* _HMP;
  int                     _ncomp;

 public:

  EdgeInfoContainer<DIM>() {}
  ~EdgeInfoContainer<DIM>();

  void BasicInit(const HierarchicalMesh*, int);
  void ModifyHanging();

  const HierarchicalMesh* GetMesh() const { return _HMP; }

  void ShowStatistics() const;
};
}

#endif
