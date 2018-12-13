/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#ifndef  __Laplace2d_h
#define  __Laplace2d_h

#include  "equation.h"
#include  "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class Laplace2d : public virtual Equation
{
  protected:
  
  mutable double visc;
  
  public:

  Laplace2d();
  Laplace2d(const ParamFile* pf);

  std::string GetName()  const { return "Laplace2d";}

  int         GetNcomp() const {return 1;}

  void OperatorStrong(DoubleVector& b, const FemFunction& U) const;
  void SetTimePattern(TimePattern& P) const;
  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};
}

#endif
