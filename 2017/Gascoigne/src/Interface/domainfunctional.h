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


#ifndef  __DomainFunctional_h
#define  __DomainFunctional_h


#include  "functional.h"
#include  "gascoigne.h"
#include  "vertex.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments DomainFunctional

  ///
  ///
  /////////////////////////////////////////////

  class DomainFunctional : public virtual Functional
  {
    private:

    protected:

    public:
      DomainFunctional() : Functional() {}
      virtual ~DomainFunctional() {}

      virtual double J(const FemFunction& U, const Vertex2d& v) const {
        std::cerr << "\"DomainFunctional::J\" not written" << std::endl; 
        abort();
      }

      virtual double J(const FemFunction& U, const Vertex3d& v) const {
        std::cerr << "\"DomainFunctional::J\" not written" << std::endl; 
        abort();
      }
  };
}

#endif
