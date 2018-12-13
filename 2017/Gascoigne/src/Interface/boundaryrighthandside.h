/**
*
* Copyright (C) 2004, 2006 by the Gascoigne 3D authors
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


#ifndef  __BoundaryRightHandSide_h
#define  __BoundaryRightHandSide_h

#include  "vertex.h"
#include  <set>
#include  "nvector.h"
#include  <string>
#include  "gascoigne.h"
#include  "application.h"

/*-------------------------------------------------------*/

namespace Gascoigne
{
  
  //////////////////////////////////////////////
  ///
  ///@brief
  /// Interface class for Boundary Conditions of Neumann
  /// or Robin type

  /// void operator()(Vector& b, const Vertex2d& v, int col)
  /// gets the coordinate v and color of boundarypart "col" and 
  /// sets the values of b. b is a vector of length ncomp
  ///
  //////////////////////////////////////////////

  class BoundaryRightHandSide : public virtual Application
  {
    private:

    protected:

    public:
      BoundaryRightHandSide() : Application() {}
      ~BoundaryRightHandSide() {}

      virtual int GetNcomp() const=0;

      virtual double operator()(int c, const Vertex2d& v, const Vertex2d& n, int color) const {
        std::cerr << "\"BoundaryRightHandSide::operator()\" not written!" << std::endl;
        abort();
      }
      virtual double operator()(int c, const Vertex3d& v, const Vertex3d& n, int color) const {
        std::cerr << "\"BoundaryRightHandSide::operator()\" not written!" << std::endl;
        abort();
      }

      virtual void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v, const Vertex2d& n, int color) const {
        for(int c=0;c<GetNcomp();c++)
        {
          b[c] += N.m()* (*this)(c,v,n,color);
        }
      }
      virtual void operator()(VectorIterator b, const TestFunction& N, const Vertex3d& v, const Vertex3d& n, int color) const {
        for(int c=0;c<GetNcomp();c++)
        {
          b[c] += N.m()* (*this)(c,v,n,color);
        }
      }

      virtual void SetCellSize(double h) const { }
  };

  typedef BoundaryRightHandSide BoundaryInitialCondition;

/*-------------------------------------------------------*/

}

#endif
