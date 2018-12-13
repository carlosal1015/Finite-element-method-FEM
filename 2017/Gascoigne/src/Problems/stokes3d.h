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


#ifndef  __Stokes3d_h
#define  __Stokes3d_h

#include  "stokes2d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class Stokes3d : public Stokes2d
{
 protected:

  double Divergence(const FemFunction& U) const;
  double Laplace(const TestFunction& U, const TestFunction& N) const;

public:

  ~Stokes3d();
  Stokes3d();
  Stokes3d(const ParamFile* pf);

  std::string GetName() const { return "Stokes3d";}

  int  GetNcomp  () const { return 4; }

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};
}

#endif
