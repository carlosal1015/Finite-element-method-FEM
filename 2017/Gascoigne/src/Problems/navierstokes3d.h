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


#ifndef  __NavierStokes3d_h
#define  __NavierStokes3d_h

#include  "navierstokes2d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class NavierStokes3d : public NavierStokes2d
{
  protected:
  
  double Laplace(const TestFunction& U, 
		 const TestFunction& N) const
    {
      return U.x()*N.x() + U.y()*N.y() + U.z()*N.z();
    }
  
  double Convection(const std::vector<TestFunction>& U, 
		    const TestFunction& N) const
    {
      return U[1].m()*N.x() + U[2].m()*N.y() + U[3].m()*N.z();
    }
  double Divergence(const std::vector<TestFunction>& U) const
    {
      return U[1].x() + U[2].y() + U[3].z();
    }
  
  public:

  ~NavierStokes3d();
  NavierStokes3d();
  NavierStokes3d(const ParamFile* pf);

  std::string GetName() const;

  int    GetNcomp  () const { return 4; }

  DoubleVector GetViscosities() const;

  void SetTimePattern(TimePattern& P) const;

  void OperatorStrong(DoubleVector& b, const FemFunction& U) const;

  void point(double h, const FemFunction& U, const Vertex3d& v) const 
    { 
      _h = h;
    }

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
  void MatrixLoop(EntryMatrix& A, const FemFunction& U, const FemFunction& M, const FemFunction& N) const;
};
}

#endif
