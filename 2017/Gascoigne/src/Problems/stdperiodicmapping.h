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


#ifndef  __stdperiodicmapping_h
#define  __stdperiodicmapping_h

#include  "periodicmapping.h"

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

  class StdPeriodicMapping : public virtual PeriodicMapping
  {
    public:
      StdPeriodicMapping() : PeriodicMapping() {}

      virtual ~StdPeriodicMapping() {}

      std::string GetName() const {return "StdPeriodicMapping";}

      virtual void transformCoords(Vertex2d& w, const Vertex2d& v) const {
      /*-------------------------------------------------------
     | Affin-lineare Transformation, die die beiden
     | periodischen Raender aufeinander abbildet.
     | v = gegebener Knoten, w = gesuchter Knoten auf dem
     | anderen Rand.
     | Die hier angegebene Identitaet funktioniert nur, wenn
     | es sich um eine Translation im rechten Winkel zum Rand
     | handelt und der Rand ein Geradenstueck ist.
     | Sonst lokal wie gewuenscht ueberschreiben.
     -------------------------------------------------------*/
        w.x() = v.x();
        w.y() = v.y();
      }
      
      virtual void transformCoords(Vertex3d& w, const Vertex3d& v) const {
      /*-------------------------------------------------------
     | Affin-lineare Transformation, die die beiden
     | periodischen Raender aufeinander abbildet.
     | v = gegebener Knoten, w = gesuchter Knoten auf dem
     | anderen Rand.
     | Die hier angegebene Identitaet funktioniert nur, wenn
     | es sich um eine Translation im rechten Winkel zum Rand
     | handelt und der Rand ein Geradenstueck ist.
     | Sonst lokal wie gewuenscht ueberschreiben.
     -------------------------------------------------------*/
        w.x() = v.x();
        w.y() = v.y();
        w.z() = v.z();
      }
      
      virtual std::set<int> preferred_colors()const {
        return std::set<int>();
      }
  };
}

#endif //__stdperiodicmapping_h

