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


#include "baseq13d.h"

/*------------------------------------------------------------------*/

namespace Gascoigne
{
BaseQ13d::BaseQ13d() 
{
  BasicInit();
}

/*------------------------------------------------------------------*/

void BaseQ13d::point(const Vertex3d& s) const
{
  for(int i=0;i<8;i++) 
    { 
      int ix = i%2;
      int iy = (i%4)/2;
      int iz = i/4;

      double px = psi(ix,s.x());
      double py = psi(iy,s.y());
      double pz = psi(iz,s.z());

      double dx = psi_x(ix,s.x());
      double dy = psi_x(iy,s.y());
      double dz = psi_x(iz,s.z());

      N  [i]     = px * py * pz;
      DN [i].x() = dx * py * pz;
      DN [i].y() = px * dy * pz;
      DN [i].z() = px * py * dz;

//       dxy[i]     = dx * dy * pz;
//       dxz[i]     = dx * py * dz;
//       dyz[i]     = px * dy * dz;
//       dxyz[i]    = dx * dy * dz;
    }
}
}

/*------------------------------------------------------------------*/

