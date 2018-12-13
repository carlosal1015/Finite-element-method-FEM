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


#ifndef  __ZeroDirichletData_h
#define  __ZeroDirichletData_h

#include  "dirichletdata.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class ZeroDirichletData : public DirichletData
{
protected:


public:

  ZeroDirichletData() : DirichletData() {}
  std::string GetName() const {return "Zero";}
  void operator()(DoubleVector& b, const Vertex2d& v, int col) const {b.zero();}
  void operator()(DoubleVector& b, const Vertex3d& v, int col) const {b.zero();}
};
}

#endif
