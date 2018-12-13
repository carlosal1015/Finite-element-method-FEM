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


#include  "q2lps2d.h"
#include  "galerkinlpsintegratorq2.h"

namespace Gascoigne
{
/* ----------------------------------------- */

void Q2Lps2d::BasicInit(const ParamFile* paramfile)
{
  if (GetIntegrator()==NULL)
    GetIntegratorPointer() =  new GalerkinLpsIntegratorQ2<2>;
  assert(GetIntegrator());

  Q22d::BasicInit(paramfile);
}

/* ----------------------------------------- */

}
