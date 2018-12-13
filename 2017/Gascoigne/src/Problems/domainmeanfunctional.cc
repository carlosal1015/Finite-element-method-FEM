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


#include  "domainmeanfunctional.h"


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
double AllDomainFunctional::J(const FemFunction& U, const Vertex2d& V) const
{
  return U[_comp].m();
}

/*-----------------------------------------*/

double AllDomainFunctional::J(const FemFunction& U, const Vertex3d& V) const
{
  return U[_comp].m();
}
/*-----------------------------------------*/

void SubDomainFunctional::SetCoordinates(double x0, double x1, double y0, double y1, double z0, double z1)
{
  _x0 = x0;
  _x1 = x1;
  _y0 = y0;
  _y1 = y1;
  _z0 = z0;
  _z1 = z1;
}

/*-----------------------------------------*/

double SubDomainFunctional::J(const FemFunction& U, const Vertex2d& V) const
{
  if ((V.x()>_x1) || (V.x()<_x0)) return 0.;
  if ((V.y()>_y1) || (V.y()<_y0)) return 0.;

  return U[_comp].m();
}

/*-----------------------------------------*/

double SubDomainFunctional::J(const FemFunction& U, const Vertex3d& V) const
{
  if ((V.x()>_x1) || (V.x()<_x0)) return 0.;
  if ((V.y()>_y1) || (V.y()<_y0)) return 0.;
  if ((V.z()>_z1) || (V.z()<_z0)) return 0.;

  return U[_comp].m();
}
}
