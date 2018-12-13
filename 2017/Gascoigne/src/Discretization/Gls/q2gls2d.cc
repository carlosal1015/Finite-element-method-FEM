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


#include  "q2gls2d.h"
#include  "galerkinglsintegratorq2.h"
#include  "transformation2d.h"
#include  "finiteelement.h"
#include  "baseq12d.h"
#include  "baseq22d.h"
#include "integratorwithsecond.h"

using namespace std;
namespace Gascoigne
{
/* ----------------------------------------- */

void Q2Gls2d::BasicInit(const ParamFile* paramfile)
{
  assert(GetIntegrator()==NULL);
  GetIntegratorPointer() =  new GalerkinGlsIntegratorQ2<2>;
  GetIntegratorPointer()->BasicInit();

  assert(GetFem()==NULL);
  typedef Transformation2d<BaseQ12d>           TransQ1;
  typedef Transformation2d<BaseQ22d>           TransQ2;
  typedef FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond> FiniteElement;
  //typedef FiniteElement<2,1,TransQ2,BaseQ22d>  FiniteElement;
  PatchDiscretization::GetFemPointer() =  new FiniteElement;

  PatchDiscretization::BasicInit(paramfile);
}
}
