/**
*
* Copyright (C) 2009, 2011 by the Gascoigne 3D authors
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


#ifndef  __periodicdata_h
#define  __periodicdata_h

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

    /// void operator()(Vertex2d& w, const Vertex2d& v, int col_w, int col_v)
    /// gets the coordinate v, color of boundarypart "col_v" and color 
    /// of boundarypart "col_w"and sets the coordinate w.
    ///
    //////////////////////////////////////////////

  class PeriodicData : public virtual Application
  {
    public:
      PeriodicData() : Application() {}

      virtual ~PeriodicData() {}

      virtual void operator()(DoubleVector& b, const Vertex2d& v, int col) const {
      /*-------------------------------------------------------
     | Falls auf einem der beiden Raender ein Offset addiert
     | werden soll, um einen Sprungterm in der Loesung zu
     | erhalten. Anwendung siehe DirichletData().
     | Nicht implementieren, falls echt periodische Loesung
     | gewuenscht ist.
     -------------------------------------------------------*/
      }

      virtual void operator()(DoubleVector& b, const Vertex3d& v, int col) const {
      /*-------------------------------------------------------
     | Falls auf einem der beiden Raender ein Offset addiert
     | werden soll, um einen Sprungterm in der Loesung zu
     | erhalten. Anwendung siehe DirichletData().
     | Nicht implementieren, falls echt periodische Loesung
     | gewuenscht ist.
     -------------------------------------------------------*/
      }
                  
      virtual std::set<int> preferred_colors()const {
        return std::set<int>();
      }
  };
}

#endif //__periodicdata_h

