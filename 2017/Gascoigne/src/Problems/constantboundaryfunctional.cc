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


#include  "constantboundaryfunctional.h"


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
ConstantBoundaryFunctional::ConstantBoundaryFunctional()
{
}

/*-----------------------------------------*/

ConstantBoundaryFunctional::~ConstantBoundaryFunctional()
{
}

/*-----------------------------------------*/

ConstantBoundaryFunctional::ConstantBoundaryFunctional
(const vector<string>& args)
{
  Construct(args);
}

/*-----------------------------------------*/

void ConstantBoundaryFunctional::Construct
(const vector<string>& args)
{
  // Baustelle
  AddColor(atoi(args[0].c_str()));
  comp  = atoi(args[1].c_str());
  value = atof(args[2].c_str());
}

/*-----------------------------------------*/

double ConstantBoundaryFunctional::J
(const FemFunction& U, const Vertex2d& v) const
{
  return value*U[comp].m();
}
}
