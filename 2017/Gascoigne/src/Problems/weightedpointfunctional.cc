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


#include "weightedpointfunctional.h"

using namespace std;

/**********************************************************/
namespace Gascoigne
{
void WeightedPointFunctional::BasicInit(const vector<Vertex2d>& v2d, const vector<int>& comps, const vector<double>& weights)
{
  _weights = weights;
  PointFunctional::BasicInit(v2d,comps);
}

/**********************************************************/

void WeightedPointFunctional::BasicInit(const vector<Vertex3d>& v3d, const vector<int>& comps, const vector<double>& weights)
{
  _weights = weights;
  PointFunctional::BasicInit(v3d,comps);
}

/**********************************************************/

double WeightedPointFunctional::J(const vector<double>& u) const
{
  double s=0;
  for(int i=0;i<u.size();++i) s+=_weights[i]*u[i];

  return s;
}
}
/**********************************************************/

