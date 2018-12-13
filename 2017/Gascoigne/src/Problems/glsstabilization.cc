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


#include  "glsstabilization.h"

/*--------------------------------------------------------------*/

namespace Gascoigne
{
GlsStabilization::GlsStabilization() : Stabilization()
{
  _delta = _tau = 0.;
  delta0 = sdelta0 = tau0 = 0.;
  _sdelta.resize(1,0.);
}

/*--------------------------------------------------------------*/

void GlsStabilization::NavierStokes(double h, double visc)
{
  _h = h;

  double val  = xeta0 * visc/(h*h) + _norm/h;
  double valc = xeta0 * visc/(h*h) + _norm/h;
  if(dt>0.)
    {
      valc += _dtfactor/dt;
    }
  _alpha = alpha0 / val;
  _delta = delta0 / valc;
  _tau   = tau0   * _norm * _norm *_delta;
}

/*--------------------------------------------------------------*/

void GlsStabilization::ConvectionDiffusion(double visc)
{
  double val = xeta0 * visc/(_h*_h) + _norm/_h;
  if(dt>0.)
    {
      val  += _dtfactor/dt;
    }
  _sdelta[0] = sdelta0 / val;
}
}
