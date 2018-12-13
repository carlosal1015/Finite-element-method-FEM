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


#include "dwrlps2d.h"

namespace Gascoigne
{

/*-------------------------------------------------*/

void DwrLps2d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
//   nmatrix<double> T;
//   for(int iq=0; iq<GetPatchMesh()->npatches(); ++iq)
//     {
//       Transformation(T,iq);
//       GetFem()->ReInit(T);
      
//       GlobalToLocal(__U,u,iq);
//       GetIntegrator()->Form(EQ,__F,*GetFem(),__U,__Q);
//       LocalToGlobal(f,__F,iq,d);
//     }
  Q1Lps2d::Form(f,u,EQ,d);
}

/*-------------------------------------------------*/

}
