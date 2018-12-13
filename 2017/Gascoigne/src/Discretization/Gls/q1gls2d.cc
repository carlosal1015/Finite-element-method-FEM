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


#include  "q1gls2d.h"
#include  "galerkinglsintegrator.h"
#include  "transformation2d.h"
#include  "finiteelement.h"
#include  "baseq12d.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
void Q1Gls2d::BasicInit(const ParamFile* pf)
{
  assert(HN==NULL);
  HN = NewHNStructure();
  assert(HN);

  assert(CellDiscretization::GetIntegratorPointer()==NULL);
  CellDiscretization::GetIntegratorPointer() =  new GalerkinGlsIntegrator<2>;

  GetIntegratorPointer()->BasicInit();

  assert(CellDiscretization::GetFemPointer()==NULL);
  typedef Transformation2d<BaseQ12d>           TransQ1;
  typedef FiniteElement<2,1,TransQ1,BaseQ12d>  FiniteElement;
  CellDiscretization::GetFemPointer() =  new FiniteElement;

  CellDiscretization::BasicInit(pf);
}
}

/* ----------------------------------------- */
