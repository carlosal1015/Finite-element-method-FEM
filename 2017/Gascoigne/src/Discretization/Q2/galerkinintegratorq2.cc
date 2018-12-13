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


#include "galerkinintegratorq2.h"

namespace Gascoigne
{

/* ----------------------------------------- */

template<int DIM>
GalerkinIntegratorQ2<DIM>::GalerkinIntegratorQ2() {}
  
/* ----------------------------------------- */

template<int DIM>
GalerkinIntegratorQ2<DIM>::~GalerkinIntegratorQ2<DIM>() {}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegratorQ2<DIM>::BasicInit()
{
  if (DIM==2)
    {
      if (!GalerkinIntegrator<DIM>::FormFormulaPointer())     GalerkinIntegrator<DIM>::FormFormulaPointer() = new QuadGauss9;
      if (!GalerkinIntegrator<DIM>::MassFormulaPointer())     GalerkinIntegrator<DIM>::MassFormulaPointer() = new QuadGauss9; // ?????
      if (!GalerkinIntegrator<DIM>::ErrorFormulaPointer())    GalerkinIntegrator<DIM>::ErrorFormulaPointer() = new QuadGauss16;
      if (!GalerkinIntegrator<DIM>::BoundaryFormulaPointer()) GalerkinIntegrator<DIM>::BoundaryFormulaPointer() = new LineGauss3;
    }
  else if (DIM==3)
    {
      if (!GalerkinIntegrator<DIM>::FormFormulaPointer())     GalerkinIntegrator<DIM>::FormFormulaPointer() = new HexGauss27;
      if (!GalerkinIntegrator<DIM>::MassFormulaPointer())     GalerkinIntegrator<DIM>::MassFormulaPointer() = new HexGauss27; // ?????
      if (!GalerkinIntegrator<DIM>::ErrorFormulaPointer())    GalerkinIntegrator<DIM>::ErrorFormulaPointer() = new HexGauss64;
      if (!GalerkinIntegrator<DIM>::BoundaryFormulaPointer()) GalerkinIntegrator<DIM>::BoundaryFormulaPointer() = new QuadGauss9;
    }
  assert(GalerkinIntegrator<DIM>::FormFormulaPointer());
  assert(GalerkinIntegrator<DIM>::ErrorFormulaPointer());
  assert(GalerkinIntegrator<DIM>::BoundaryFormulaPointer());
  assert(GalerkinIntegrator<DIM>::MassFormulaPointer());
}

/* ----------------------------------------- */

template class GalerkinIntegratorQ2<2>;
template class GalerkinIntegratorQ2<3>;

}
