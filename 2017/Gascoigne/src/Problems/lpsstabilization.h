/**
*
* Copyright (C) 2004, 2006 by the Gascoigne 3D authors
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


#ifndef __LpsStabilization_h
#define __LpsStabilization_h

#include  "stabilization.h"

/*-------------------------------------------*/

namespace Gascoigne
{
class LpsStabilization : public Stabilization
{
 protected:

  double _delta, _tau;
  nvector<double> _sdelta;

  void NavierStokes(double h, double visc);

 public:

  LpsStabilization();

  double delta0, sdelta0, tau0;

  double&  tau()      { return _tau;}
  double   tau()      const { return _tau;}
  double&  delta()      { return _delta;}
  double  delta()      const { return _delta;}
  double  delta(int i) const { return _sdelta[i];}
  void BasicInit(int n);
  void ReInit(double h, double visc, double u, double v)
    {
      norm(u,v);
      NavierStokes(h,visc);
    };
  void ReInit(double h, double visc, double u, double v, double w)
    {
      norm(u,v,w);
      NavierStokes(h,visc);
    };
  void ConvectionDiffusion(double visc);
  void ConvectionDiffusion(const nvector<double>& visc);
};
}

/*-------------------------------------------*/

#endif
