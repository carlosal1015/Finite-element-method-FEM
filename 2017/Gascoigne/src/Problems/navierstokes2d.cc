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


#include  "navierstokes2d.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{


NavierStokes2d::~NavierStokes2d()
{
}

/*-----------------------------------------*/

NavierStokes2d::NavierStokes2d() : Equation()
{
  _penalty = 0.; _visc = 0.01; _h = 0.;
}

/*-----------------------------------------*/

NavierStokes2d::NavierStokes2d(const ParamFile* pf) : Equation()
{
  _h = 0.;
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 1.);
  DFH.insert("cut"  , &_cut  ,  1.e10);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

void NavierStokes2d::OperatorStrong(DoubleVector& b, const FemFunction& U) const
{
  b[0] = Divergence(U);
  b[1] = Convection(U,U[1]) - _visc * U[1].D() + U[0].x();
  b[2] = Convection(U,U[2]) - _visc * U[2].D() + U[0].y();
}

/*-----------------------------------------*/

double NavierStokes2d::Laplace(const TestFunction& U, const TestFunction& N) const
{
  return U.x()*N.x() + U.y()*N.y();
}

/*-----------------------------------------*/

double NavierStokes2d::Convection(const FemFunction& U, const TestFunction& N) const
{
  return U[1].m()*N.x() + U[2].m()*N.y();
}

/*-----------------------------------------*/

double NavierStokes2d::Divergence(const FemFunction& U) const
{
  return U[1].x() + U[2].y();
}
 
/*-----------------------------------------*/

void NavierStokes2d::SetTimePattern(TimePattern& P) const
{
  P(0,0) = _penalty;
  P(1,1) = 1.;
  P(2,2) = 1.;
}

/*-----------------------------------------*/

void NavierStokes2d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  ////////////// Continuity ////////////////////////////////////////////////

  b[0] += Divergence(U) * N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  b[1] -= U[0].m()*N.x();
  b[2] -= U[0].m()*N.y();
	  
  b[1] += Convection(U,U[1]) * N.m();
  b[2] += Convection(U,U[2]) * N.m();

  // viscous terms
  b[1] += _visc * Laplace(U[1],N);
  b[2] += _visc * Laplace(U[2],N);
}

/*-----------------------------------------*/

void NavierStokes2d::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  double MN = M.m()*N.m();
  double Mx = M.x()*N.m();
  double My = M.y()*N.m();
  double laplace = Laplace(M,N);
  
  ////////////// Continuity ////////////////////////////////////////////////

  A(0,1) += Mx;
  A(0,2) += My;
  
  ////////////// Momentum ////////////////////////////////////////////////

  A(1,0) -= M.m()*N.x();
  A(2,0) -= M.m()*N.y();
  
  double cl = Convection(U,M) * N.m() + _visc*laplace;
  
  A(1,1) += cl;
  A(2,2) += cl;
  
  double cut = _cut * _h;
  
  A(1,1) += Gascoigne::max(U[1].x()*MN, -cut);
  A(2,2) += Gascoigne::max(U[2].y()*MN, -cut);
  A(1,2) += Gascoigne::min(U[1].y()*MN, cut);
  A(2,1) += Gascoigne::min(U[2].x()*MN, cut);
}
}
/*-----------------------------------------*/
