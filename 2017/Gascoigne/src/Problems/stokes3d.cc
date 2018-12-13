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


#include  "stokes3d.h"

#define P U[0]

/*-----------------------------------------*/

namespace Gascoigne
{
Stokes3d::~Stokes3d()
{
}

/*-----------------------------------------*/

Stokes3d::Stokes3d() : Stokes2d()
{}

/*-----------------------------------------*/

Stokes3d::Stokes3d(const ParamFile* pf) : Stokes2d(pf)
{
}

/*-----------------------------------------*/

double Stokes3d::Divergence(const FemFunction& U) const
{
  return U[1].x() + U[2].y() + U[3].z();
}

/*-----------------------------------------*/

double Stokes3d::Laplace(const TestFunction& U, const TestFunction& N) const
{
  return U.x()*N.x() + U.y()*N.y() + U.z()*N.z();
}

/*-----------------------------------------*/

void Stokes3d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  ////////////// Continuity //////////////////////////

  b[0] += Divergence(U) * N.m();

  ////////////// Momentum ////////////////////////////

  b[1] -= P.m()*N.x();
  b[2] -= P.m()*N.y();
  b[3] -= P.m()*N.z();
	  
  // viscous terms
  b[1] += _visc * Laplace(U[1],N);
  b[2] += _visc * Laplace(U[2],N);
  b[3] += _visc * Laplace(U[3],N);
}

/*-----------------------------------------*/

void Stokes3d::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  double laplace = Laplace(M,N);
     
  ////////////// Continuity /////////////////////////////

  A(0,0) += M.m()*N.m() * _penalty;
  A(0,1) += M.x()*N.m();
  A(0,2) += M.y()*N.m();
  A(0,3) += M.z()*N.m();

  ////////////// Momentum ///////////////////////////////

  A(1,0) -= M.m()*N.x();
  A(2,0) -= M.m()*N.y();
  A(3,0) -= M.m()*N.z();

  A(1,1) += _visc*laplace;
  A(2,2) += _visc*laplace;
  A(3,3) += _visc*laplace;
}
}

/*-----------------------------------------*/

#undef P
