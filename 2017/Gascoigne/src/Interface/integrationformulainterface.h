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


#ifndef __integrationformulainterface_h
#define __integrationformulainterface_h

#include  "nvector.h"
#include  "vertex.h"

/*------------------------------------------------------------*/

namespace Gascoigne
{
  class IntegrationFormulaInterface
  {
    public:
      IntegrationFormulaInterface() {}
      virtual ~IntegrationFormulaInterface() {}

      virtual int    n()      const=0;
      virtual double w(int k) const=0;

      virtual void xi(Vertex1d& v, int k) const {
        std::cerr << "\"IntegrationFormulaInterface::xi\" not written!" << std::endl;
        abort();
      }
      virtual void xi(Vertex2d& v, int k) const {
        std::cerr << "\"IntegrationFormulaInterface::xi\" not written!" << std::endl;
        abort();
      }
      virtual void xi(Vertex3d& v, int k) const {
        std::cerr << "\"IntegrationFormulaInterface::xi\" not written!" << std::endl;
        abort();
      }
  };
}

/*------------------------------------------------------------*/

#endif
