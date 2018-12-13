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


#include  "transformation2d.h"
#include  "baseq22dwithsecond.h"
#include  "../Q1/finiteelement.xx"
#include  "finiteelementwithsecond.xx"
#include  "finiteelementwithsecond2d.xx"

/*-----------------------------------------------------*/

namespace Gascoigne
{
template class FiniteElement          <2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond>;
template class FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond>;
}
/*-----------------------------------------------------*/

