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


#ifndef __GlsStabilization_h
#define __GlsStabilization_h

#include  "stabilization.h"
#include  "gascoigne.h"

/*-------------------------------------------*/

namespace Gascoigne
{
class GlsStabilization : public Stabilization
{
 protected:

  double _delta, _tau;
  DoubleVector _sdelta;

  void NavierStokes(double h, double visc);

 public:

  GlsStabilization();

  double delta0, tau0, sdelta0;

  double  delta()      const { return _delta;}
  double& delta()            { return _delta;}
  double  delta(int i) const { return _sdelta[i];}
  double  tau()        const { return _tau;}
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
};
}

#endif
