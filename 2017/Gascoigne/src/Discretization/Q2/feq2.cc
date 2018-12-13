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


#include  "finiteelement.h"
#include  "../Q1/finiteelement.xx"
#include  "transformation2d.h"
#include  "transformation3d.h"
#include  "baseq12d.h"
#include  "baseq22d.h"
#include  "baseq23d.h"
#include  "baseq1patch.h"
#include  "baseq13dpatch.h"

namespace Gascoigne
{
/*-----------------------------------------------------*/

typedef Transformation2d<BaseQ12d>  TQ1_2D;
typedef Transformation2d<BaseQ22d>  TQ2_2D;

template class FiniteElement<2,1,TQ1_2D,BaseQ22d>;
template class FiniteElement<2,1,TQ2_2D,BaseQ22d>;
template class FiniteElement<2,1,TQ1_2D,BaseQ12dPatch>;

/*-----------------------------------------------------*/

typedef Transformation3d<BaseQ13d>  TQ1_3D;
typedef Transformation3d<BaseQ23d>  TQ2_3D;

template class FiniteElement<3,2,TQ2_3D,BaseQ23d>;
template class FiniteElement<3,2,TQ1_3D,BaseQ13dPatch>;
}
/*-----------------------------------------------------*/
