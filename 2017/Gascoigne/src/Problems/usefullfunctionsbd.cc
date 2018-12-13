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


#include "usefullfunctionsbd.h"
#include <math.h>

/**********************************************************/

namespace Gascoigne
{
double ParabelFunction(double x, double n0, double n1)
{
  return -4.*(x-n0)*(x-n1)/((n0-n1)*(n0-n1));
}

double GaussFunction(double x, double a, double stiff)
{
  return exp( -stiff*(x-a)*(x-a) );
}

double StiffnessFunction(double x, double a, double stiff)
{
  return 0.5 * (1.+tanh(stiff*(x-a)));
}

double PlugFlowFunction(double x, double a, double b, double stiff)
{
  double h = fabs(a-b);

  return tanh(stiff/h*(x-a)) * tanh(stiff/h*(-x+b));
}
}

/**********************************************************/

