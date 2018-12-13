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


#include  "navierstokesgls2d.h"
#include  "filescanner.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
NavierStokesGls2d::~NavierStokesGls2d()
{
}

/*-----------------------------------------*/

NavierStokesGls2d::NavierStokesGls2d() 
  : GlsEquation(), NavierStokes2d() 
{
  _penalty = 0.; _visc = 0.01;
  
  ST.delta0 = ST.alpha0 = 0.2;
  ST.xeta0 = 6.; 
}

/*-----------------------------------------*/

NavierStokesGls2d::NavierStokesGls2d(const ParamFile* pf) 
  : GlsEquation(), NavierStokes2d() 
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("cut"  , &_cut  ,  1.e8);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("delta", &ST.delta0, 0.25);
  DFH.insert("tau"  , &ST.tau0  , 0.);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

void NavierStokesGls2d::glspoint(double h, const FemFunction& U, const Vertex2d& v)const
{
  ST.ReInit(h,_visc,U[1].m(),U[2].m());
}

/*-----------------------------------------*/

void NavierStokesGls2d::L(DoubleVector& dst, const FemFunction& U) const
{
  dst[0] = Divergence(U);
  dst[1] = Convection(U,U[1]) + U[0].x();
  dst[2] = Convection(U,U[2]) + U[0].y();
}

/*-----------------------------------------*/

void NavierStokesGls2d::S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const
{
  dst(0,1) = ST.alpha() * N.x();
  dst(0,2) = ST.alpha() * N.y();
  dst(1,0) = ST.tau()   * N.x();
  dst(2,0) = ST.tau()   * N.y();

  dst(1,1) = ST.delta() * Convection(U,N);
  dst(2,2) = ST.delta() * Convection(U,N);
}

/*-----------------------------------------*/

void NavierStokesGls2d::LMatrix(nmatrix<double>& A, const FemFunction& U, const TestFunction& N) const
{
  A(0,1) = N.x();
  A(0,2) = N.y();
     
  A(1,0) = N.x();
  A(2,0) = N.y();

  double cl = Convection(U,N);

  A(1,1) = cl + U[1].x()*N.m();
  A(2,2) = cl + U[2].y()*N.m();
  A(1,2) =      U[1].y()*N.m();
  A(2,1) =      U[2].x()*N.m();
}

/*-----------------------------------------*/

void NavierStokesGls2d::SMatrix(DoubleVector& dst, const FemFunction& U, const FemFunction& M, const FemFunction& N) const
{
  dst[1] = ST.delta() * Convection(M,N[1]);
  dst[2] = ST.delta() * Convection(M,N[2]);
}
}
