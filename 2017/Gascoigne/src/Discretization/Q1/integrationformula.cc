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


#include  "integrationformula.h"
#include  "integrationformulasummed.h"
#include  "tensorformula2d.xx"
#include  "tensorformula3d.xx"

/*------------------------------------------------------------*/

namespace Gascoigne
{
template class IntegrationFormulaBase<1>;
template class IntegrationFormulaBase<2>;
template class IntegrationFormulaBase<3>;

/*------------------------------------------------------------*/

template class TensorFormula2d<1,LineMidPoint>;
template class TensorFormula2d<2,LineTrapez  >;
template class TensorFormula2d<3,LineSimpson >;
template class TensorFormula2d<1,LineGauss1  >;
template class TensorFormula2d<2,LineGauss2  >;
template class TensorFormula2d<3,LineGauss3  >;
template class TensorFormula2d<4,LineGauss4  >;
template class TensorFormula2d<5,LineGauss5  >;
template class TensorFormula2d<6,LineGauss6  >;
template class TensorFormula2d<7,LineGauss7  >;
template class TensorFormula2d<8,LineGauss8  >;
template class TensorFormula2d<9,LineGauss9  >;
template class TensorFormula2d<10,LineGauss10>;

template class TensorFormula3d<1,LineGauss1>;
template class TensorFormula3d<2,LineTrapez>;
template class TensorFormula3d<2,LineGauss2>;
template class TensorFormula3d<3,LineGauss3>;
template class TensorFormula3d<4,LineGauss4>;
template class TensorFormula3d<5,LineGauss5>;

/*------------------------------------------------------------*/
template class IntegrationFormulaSummed2d<QuadGauss1>;
template class IntegrationFormulaSummed2d<QuadGauss4>;
template class IntegrationFormulaSummed3d<QuadGauss1>;
template class IntegrationFormulaSummed2d<QuadTrapez>;

template class IntegrationFormulaSummed3d<HexGauss8>;
template class IntegrationFormulaSummed3d<HexTrapez>;
}

