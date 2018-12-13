/**
*
* Copyright (C) 2009 by the Gascoigne 3D authors
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


#ifndef  __periodicmapping_h
#define  __periodicmapping_h

#include  "vertex.h"
#include  <set>
#include  "nvector.h"
#include  <string>
#include  "application.h"


/*-----------------------------------------*/

namespace Gascoigne
{

    //////////////////////////////////////////////
    ///
    ///@brief
    /// Periodic Boundary Conditions

    /// void operator()(Vertex2d& w, const Vertex2d& v)
    /// gets the coordinate v of a vertex and sets the
    /// coordinate w of corresponding vertex on other
    /// boundary.
    ///
    //////////////////////////////////////////////

  class PeriodicMapping : public virtual Application
  {
    public:
      PeriodicMapping() : Application() {}

      virtual ~PeriodicMapping() {}

      virtual void transformCoords(Vertex2d& w, const Vertex2d& v) const = 0;
      virtual void transformCoords(Vertex3d& w, const Vertex3d& v) const = 0;

      virtual std::set<int> preferred_colors()const {
        return std::set<int>();
      }
  };
}

#endif //__periodicmapping_h

