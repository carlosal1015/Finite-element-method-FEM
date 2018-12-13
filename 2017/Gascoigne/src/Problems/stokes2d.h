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


#ifndef  __Stokes2d_h
#define  __Stokes2d_h

#include  "equation.h"
#include  "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class Stokes2d : public virtual Equation
{
protected:

  double _visc;
  double _penalty;

  double Laplace(const TestFunction& U, const TestFunction& N) const;
  double Divergence(const FemFunction& U) const;

public:

  ~Stokes2d();
  Stokes2d();
  Stokes2d(const ParamFile* pf);

  std::string GetName() const { return "Stokes2d";}

  int GetNcomp  () const { return 3; }

  //
  // Time
  //

  void SetTimePattern(TimePattern& P) const;

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};
}

#endif
