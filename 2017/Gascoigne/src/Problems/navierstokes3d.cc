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


#include  "navierstokes3d.h"


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
NavierStokes3d::~NavierStokes3d() {}

/*-----------------------------------------*/

NavierStokes3d::NavierStokes3d() : NavierStokes2d() {}

/*-----------------------------------------*/

void NavierStokes3d::SetTimePattern(TimePattern& P) const
{
  P(0,0) = _penalty;
  P(1,1) = 1.;
  P(2,2) = 1.;
  P(3,3) = 1.;
}

/*-----------------------------------------*/

NavierStokes3d::NavierStokes3d(const ParamFile* pf) 
  : NavierStokes2d(pf) {}


/*-----------------------------------------*/

void NavierStokes3d::OperatorStrong(DoubleVector& b, const FemFunction& U) const
{
  b[0] = Divergence(U);
  b[1] = Convection(U,U[1]) - _visc * U[1].D() + U[0].x();
  b[2] = Convection(U,U[2]) - _visc * U[2].D() + U[0].y();
  b[3] = Convection(U,U[3]) - _visc * U[3].D() + U[0].z();
}

/*-----------------------------------------*/

void NavierStokes3d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  ////////////// Continuity ////////////////////////////////////////////////

  b[0] += Divergence(U) * N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  b[1] -= U[0].m()*N.x();
  b[2] -= U[0].m()*N.y();
  b[3] -= U[0].m()*N.z();
	  
  b[1] += Convection(U,U[1]) * N.m();
  b[2] += Convection(U,U[2]) * N.m();
  b[3] += Convection(U,U[3]) * N.m();

  // viscous terms
  b[1] += _visc * Laplace(U[1],N);
  b[2] += _visc * Laplace(U[2],N);
  b[3] += _visc * Laplace(U[3],N);
}

/*-----------------------------------------*/

void NavierStokes3d::Matrix
(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  double MN = M.m()*N.m();
     
  ////////////// Continuity ////////////////////////////////////////////////

  A(0,0) += MN * _penalty;
  A(0,1) += M.x()*N.m();
  A(0,2) += M.y()*N.m();
  A(0,3) += M.z()*N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  A(1,0) -= M.m()*N.x();
  A(2,0) -= M.m()*N.y();
  A(3,0) -= M.m()*N.z();

  double sum = Convection(U,M)*N.m() + _visc*Laplace(M,N);
  A(1,1) += sum;
  A(2,2) += sum;
  A(3,3) += sum;

  double tau = _cut * _h;

  A(1,1) += Gascoigne::max(U[1].x()*MN, -tau);
  A(2,2) += Gascoigne::max(U[2].y()*MN, -tau);
  A(3,3) += Gascoigne::max(U[3].z()*MN, -tau);

  A(1,2) += Gascoigne::min(U[1].y()*MN, tau);
  A(1,3) += Gascoigne::min(U[1].z()*MN, tau);
  A(2,1) += Gascoigne::min(U[2].x()*MN, tau);
  A(2,3) += Gascoigne::min(U[2].z()*MN, tau);
  A(3,1) += Gascoigne::min(U[3].x()*MN, tau);
  A(3,2) += Gascoigne::min(U[3].y()*MN, tau);
}

/*-----------------------------------------*/

string NavierStokes3d::GetName() const 
{ 
  return "NavierStokes3d";
}
}

/*-----------------------------------------*/
