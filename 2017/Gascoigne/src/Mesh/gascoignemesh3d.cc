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


#include  "gascoignemesh3d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
GascoigneMesh3d::GascoigneMesh3d()
{
}

/*---------------------------------------------------*/

IntVector GascoigneMesh3d::IndicesOfCell(int iq) const
{
  IntVector indices(8);
  
  int offset = 8*iq;

  indices[0] = nc[offset];
  indices[1] = nc[offset+1];
  indices[2] = nc[offset+3];
  indices[3] = nc[offset+2];
  indices[4] = nc[offset+4];
  indices[5] = nc[offset+5];
  indices[6] = nc[offset+7];
  indices[7] = nc[offset+6];

  return indices;
}
}
