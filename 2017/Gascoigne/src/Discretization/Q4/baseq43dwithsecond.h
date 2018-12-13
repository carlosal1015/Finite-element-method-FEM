/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#ifndef __baseq43dwithsecond_h
#define __baseq43dwithsecond_h

#include  "baseq43d.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  computes in addition second derivatives

///
///
/////////////////////////////////////////////

class BaseQ43dWithSecond : public BaseQ43d
{

 protected:

  mutable DoubleVector dxx,dyy,dzz,dxy,dxz,dyz;

 public:

  BaseQ43dWithSecond();

  void point(const Vertex3d& s) const;

  double phi_xx(int i) const {return dxx[i];}
  double phi_yy(int i) const {return dyy[i];}
  double phi_zz(int i) const {return dzz[i];}
  double phi_xy(int i) const {return dxy[i];}
  double phi_xz(int i) const {return dxz[i];}
  double phi_yz(int i) const {return dyz[i];}
};
}

#endif
