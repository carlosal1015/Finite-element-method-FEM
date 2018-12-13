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


#ifndef  __GascoigneMeshTransferConstructor_h
#define  __GascoigneMeshTransferConstructor_h

#include  "gascoignemeshtransfer.h"
#include  "levelmesh2d.h"
#include  "levelmesh3d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class GascoigneMeshTransferConstructor2d
{
public:

  GascoigneMeshTransferConstructor2d(const HierarchicalMesh2d* HM, GascoigneMeshTransfer* GMT, 
				     const LevelMesh2d* LMfine, const LevelMesh2d* LMcoarse);
};

/*-----------------------------------------*/

class GascoigneMeshTransferConstructor3d
{
public:

  GascoigneMeshTransferConstructor3d(const HierarchicalMesh3d* HM, GascoigneMeshTransfer* GMT, 
				     const LevelMesh3d* LMfine, const LevelMesh3d* LMcoarse);
};
}

/*-----------------------------------------*/

#endif
