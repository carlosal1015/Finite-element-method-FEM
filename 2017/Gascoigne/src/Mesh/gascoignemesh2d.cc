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


#include  "gascoignemesh2d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
GascoigneMesh2d::GascoigneMesh2d()
{
}

/*-----------------------------------------*/

GascoigneMesh2d::~GascoigneMesh2d() 
{
}

/*---------------------------------------------------*/

IntVector GascoigneMesh2d::IndicesOfCell(int iq) const
{
  IntVector indices(4);

  int iq4 = iq*4;

  indices[0] = nc[iq4];
  indices[1] = nc[iq4+1];
  indices[2] = nc[iq4+3];
  indices[3] = nc[iq4+2];

  return indices;
}
}
