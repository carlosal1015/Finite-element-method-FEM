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


#ifndef  __BechmarkFunctionals_h
#define  __BechmarkFunctionals_h

#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"

/*-----------------------------------------*/
/* Functionals for the 2d Navier-Stokes Benchmark*/
/*-----------------------------------------*/

namespace Gascoigne
{

class DragFunctional : public virtual Gascoigne::ResidualFunctional
{
  public:

  DragFunctional() : ResidualFunctional() 
    {
    __comps.push_back(1);
    __cols.insert(80);
    __scales.push_back(50);
    ExactValue() = 5.579535;       // fuer den runden
;

    __DD  = new Gascoigne::DirichletDataByColor(GetComps(),GetColors(),GetScales());
  }

  std::string GetName() const { return "DragFunctional";}
};

}

/*-----------------------------------------*/

#endif
