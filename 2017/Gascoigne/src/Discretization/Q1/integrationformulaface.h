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


#ifndef  __IntegrationFormulaFace_h
#define  __IntegrationFormulaFace_h

#include  "integrationformula.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class IntegrationFormulaEdge1 : public IntegrationFormula2d
{
protected:
 public:  
  IntegrationFormulaEdge1() : IntegrationFormula2d(4)
    {
      iw = 0.25;
      ic[0].x() = 0.5;  ic[0].y() = 0.0;
      ic[1].x() = 1.0;  ic[1].y() = 0.5;
      ic[2].x() = 0.5;  ic[2].y() = 1.0;
      ic[3].x() = 0.0;  ic[3].y() = 0.5;
    }

};
}

#endif
