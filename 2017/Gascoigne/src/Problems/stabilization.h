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


#ifndef __stabilization_h
#define __stabilization_h

#include "nvector.h"

/*-------------------------------------------*/

namespace Gascoigne
{
  class Stabilization
{
 protected:

  double _alpha, _h, dt, _dtfactor, _norm;

  void norm(double u, double v)
    { 
      _norm = sqrt(u*u+v*v)+1e-6; 
    }
  void norm(double u, double v, double w) 
    { 
      _norm = sqrt(u*u+v*v+w*w)+1e-6; 
    }

 public:

  Stabilization() : dt(0.), _dtfactor(1.), _norm(0.), alpha0(0.2), xeta0(6.) {}
  ~Stabilization() {}

  double alpha0, xeta0;

  double& DeltaT()           { return dt;}
  double  DeltaT()     const { return dt;}
  double  alpha()      const { return _alpha;}
  double& alpha()            { return _alpha;}
  double  h()          const { return _h;}
  double& dtfactor()         { return _dtfactor; }
  void    SetConvection(double u, double v)           { norm(u,v);}    
  void    SetConvection(double u, double v, double w) { norm(u,v,w);}    

  void ReInit(double h, double visc)
    {
      _h = h;
      _alpha = alpha0 *h*h / (xeta0 * visc);
    }
};
}

#endif
