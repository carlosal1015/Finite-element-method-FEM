/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#include  "baseq22dwithsecond.h"


using namespace std;

namespace Gascoigne{

#define NDOF   9
#define NDOF1d 3

/*-----------------------------------------*/

BaseQ22dWithSecond::BaseQ22dWithSecond() : BaseQ22d()
{
  dxx.reservesize(NDOF);
  dxy.reservesize(NDOF);
  dyy.reservesize(NDOF);
}

/*------------------------------------------------------------------*/

void BaseQ22dWithSecond::point(const Vertex2d& s) const 
{
  BaseQ22d::point(s);
  for(int i=0;i<NDOF;i++) 
    { 
      int ix = i%NDOF1d;
      int iy = i/NDOF1d;

      dxx[i]     = psi_xx(ix,s.x()) * psi   (iy,s.y());
      dxy[i]     = psi_x (ix,s.x()) * psi_x (iy,s.y());
      dyy[i]     = psi   (ix,s.x()) * psi_xx(iy,s.y());
    }
}

/*------------------------------------------------------------------*/

#undef NDOF
#undef NDOF1d

}
