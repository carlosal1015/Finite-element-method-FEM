/**
*
* Copyright (C) 2004, 2005, 2006, 2007 by the Gascoigne 3D authors
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


#include  "patchintegrationformula.h"
#include  "patchintegrationformula.xx"

namespace Gascoigne
{

/*------------------------------------------------------------*/

template class PatchFormula1d<2,LineGauss2>;
template class PatchFormula1d<3,LineGauss3>;
template class PatchFormula1d<4,LineGauss4>;

/*------------------------------------------------------------*/

template class PatchFormula2d<4,TensorFormula2d<2,LineTrapez> >;
template class PatchFormula2d<16,TensorFormula2d<4,LineGauss4> >;
template class PatchFormula2d<9,TensorFormula2d<3,LineGauss3> >;
template class PatchFormula2d<4,TensorFormula2d<2,LineGauss2> >;

/*------------------------------------------------------------*/

template class PatchFormula3d<8, TensorFormula3d<2,LineGauss2> >;
template class PatchFormula3d<27,TensorFormula3d<3,LineGauss3> >;
template class PatchFormula3d<64,TensorFormula3d<4,LineGauss4> >;
template class PatchFormula3d<125,TensorFormula3d<5,LineGauss5> >;

/*------------------------------------------------------------*/

}

