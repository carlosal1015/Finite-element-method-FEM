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


#include  "laplace3d.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{
Laplace3d::Laplace3d(const ParamFile* pf) : 
  Laplace2d(pf)
{
  DataFormatHandler DFH;
  DFH.insert("betax",&betax,0.);
  DFH.insert("betay",&betay,0.);
  DFH.insert("betaz",&betaz,0.);
  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

void Laplace3d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  b[0] += visc* (U[0].x()*N.x()+U[0].y()*N.y()+U[0].z()*N.z());
  b[0] += (betax * U[0].x() + betay * U[0].y() + betaz * U[0].z()) * N.m();
}

/*-----------------------------------------*/

void Laplace3d::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  A(0,0) += visc* (M.x()*N.x()+M.y()*N.y()+M.z()*N.z());
  A(0,0) += (betax * M.x() + betay * M.y() + betaz * M.z()) * N.m();
}


/*-----------------------------------------*/

void Laplace3d::OperatorStrong(DoubleVector& b, const FemFunction& U)const
{
  b[0] -= visc*U[0].D();
  b[0] += betax * U[0].x() + betay * U[0].y() + betaz * U[0].z();
}
}
