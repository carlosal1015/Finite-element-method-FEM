/**
*
* Copyright (C) 2008 by the Gascoigne 3D authors
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


#include "dynamicblockmatrix.xx"
#include "fmatrixblock.h"
#include "cfdblock3d.h"

namespace Gascoigne
{
template class DynamicBlockMatrix<FMatrixBlock<1> >;
template class DynamicBlockMatrix<FMatrixBlock<2> >;
template class DynamicBlockMatrix<FMatrixBlock<3> >;
template class DynamicBlockMatrix<FMatrixBlock<4> >;


template class DynamicBlockMatrix<CFDBlock3d>;
}
