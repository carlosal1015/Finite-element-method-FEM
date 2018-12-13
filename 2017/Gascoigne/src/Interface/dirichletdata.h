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


#ifndef  __dirichletdata_h
#define  __dirichletdata_h

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
  /// Interface class for Dirichlet Boundary Conditions

  /// void operator()(Vector& b, const Vertex2d& v, int col)
  /// gets the coordinate v and color of boundarypart "col" and 
  /// sets the values of b. b is a vector of length ncomp
  ///
  //////////////////////////////////////////////

  class DirichletData : public virtual Application
  {
    private:

    protected:

    public:
      DirichletData() : Application() {}
      virtual ~DirichletData() {}

      virtual void operator()(DoubleVector& b, const Vertex2d& v, int col) const {
        std::cerr << "\"DirichletData::operator()\" not written!" << std::endl;
        abort();
      }

      virtual void operator()(DoubleVector& b, const Vertex3d& v, int col) const {
        std::cerr << "\"DirichletData::operator()\" not written!" << std::endl;
        abort();
      }

      virtual std::set<int> preferred_colors()const {
        return std::set<int>();
      }
  };
}

#endif
