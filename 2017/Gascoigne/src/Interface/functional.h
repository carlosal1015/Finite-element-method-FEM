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


#ifndef  __Functional_h
#define  __Functional_h

#include  "application.h"


/*-----------------------------------------*/

namespace Gascoigne
{

  //////////////////////////////////////////////
  //
  /// Interface class for Functional
  //
  //////////////////////////////////////////////

  class Functional : public virtual Application
  {
    private:
      
    protected:
      double  exact;
      bool    exactisknown;

    public:
      Functional() : Application(), exact(0.), exactisknown(0) {}
      ~Functional() {}
      Functional(const Functional& F) : Application(F) {
        exact = F.ExactValue();
      } 

      double  ExactValue() const { 
        return exact;
      }
      double& ExactValue() {
        exactisknown = 1; 
        return exact;
      }
      bool ExactValueIsKnown() const 
      { 
        return exactisknown; 
      }
  };
}

#endif
