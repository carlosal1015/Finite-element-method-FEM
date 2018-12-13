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


#ifndef __baseq42dwithsecond_h
#define __baseq42dwithsecond_h

#include  "baseq42d.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  computes in addition second derivatives

///
///
/////////////////////////////////////////////

/**************************************************/

class BaseQ42dWithSecond : public BaseQ42d
{
 protected:

  mutable DoubleVector dxx, dxy, dyy;

 public:

  BaseQ42dWithSecond();

  void point(const Vertex2d& s) const;

  double phi_xx(int i) const { return dxx[i];}
  double phi_yy(int i) const { return dyy[i];}
  double phi_xy(int i) const { return dxy[i];}
};
}

#endif
