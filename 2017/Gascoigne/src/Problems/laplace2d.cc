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


#include  "laplace2d.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{
Laplace2d:: Laplace2d()
{
  visc = 1.;
}

/*-----------------------------------------*/

Laplace2d::Laplace2d(const ParamFile* pf) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("visc",&visc,1.);
  FileScanner FS(DFH,pf,"Equation");
}
 
/*-----------------------------------------*/

void Laplace2d::OperatorStrong(DoubleVector& b, const FemFunction& U)const
{
  b[0] -= visc*U[0].D();
}
 
/*-----------------------------------------*/

void Laplace2d::SetTimePattern(TimePattern& TP) const
{
  TP(0,0) = 1.;
}

/*-----------------------------------------*/

void Laplace2d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  b[0] += visc* (U[0].x()*N.x()+U[0].y()*N.y());
}

/*-----------------------------------------*/

void Laplace2d::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  A(0,0) += visc* (M.x()*N.x()+M.y()*N.y());
}
}

/*-----------------------------------------*/
