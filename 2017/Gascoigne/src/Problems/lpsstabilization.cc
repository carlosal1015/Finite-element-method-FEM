/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#include  "lpsstabilization.h"

/*--------------------------------------------------------------*/

namespace Gascoigne
{
LpsStabilization::LpsStabilization() : Stabilization(), _sdelta(0)
{
  _delta = delta0 = sdelta0 = tau0 =  0.;
}

/*--------------------------------------------------------------*/

void LpsStabilization::BasicInit(int n)
{
  _sdelta.resize(n,0.);
}

/*--------------------------------------------------------------*/

void LpsStabilization::NavierStokes(double h, double visc)
{
  _h = h;

  double val  = xeta0 * visc/(h*h) + _norm/h;
  double valc = xeta0 * visc/(h*h) + _norm/h;

  if(dt>0.)
    {
      //val  += _dtfactor/dt;
      valc += _dtfactor/dt;
    }
  _alpha = alpha0 / val;
  _delta = delta0 / valc;
  _tau   = tau0 / valc;
}

/*--------------------------------------------------------------*/

void LpsStabilization::ConvectionDiffusion(double visc)
{
  double val = xeta0 * visc/(_h*_h) + _norm/_h;
  if(dt>0.)
    {
      val  += _dtfactor/dt;
    }
  _sdelta[0] = sdelta0 / val; 
}

/*--------------------------------------------------------------*/

void LpsStabilization::ConvectionDiffusion(const nvector<double>& visc)
{
  double val = _norm/_h;
  if(dt>0.)
    {
      val += _dtfactor/dt;
    }
  for (int i=0; i<_sdelta.size(); i++)
    {
      _sdelta[i] = sdelta0 / (val+xeta0 * visc[i]/(_h*_h) ); 
    }
}
}
