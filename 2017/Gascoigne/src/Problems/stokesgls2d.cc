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


#include  "stokesgls2d.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{
StokesGls2d::~StokesGls2d()
{
}

/*-----------------------------------------*/

StokesGls2d::StokesGls2d() : GlsEquation(), Stokes2d() 
{
  _penalty = 0.; 
  _visc = 1.;
  ST.alpha0 = 1.;
}
 
/*-----------------------------------------*/

StokesGls2d::StokesGls2d(const ParamFile* pf) : GlsEquation(), Stokes2d() 
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 1.);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

void StokesGls2d::glspoint(double h, const FemFunction& U, const Vertex2d& v)const
{
  ST.ReInit(h,_visc);
}

/*-----------------------------------------*/

void StokesGls2d::L(DoubleVector& dst, const FemFunction& U) const
{
  dst[0] = Divergence(U);
  dst[1] = U[0].x();
  dst[2] = U[0].y();
}

/*-----------------------------------------*/

void StokesGls2d::S(nmatrix<double>& dst, const FemFunction& U, 
		  const TestFunction& N) const
{
  dst(0,1) = ST.alpha() * N.x();
  dst(0,2) = ST.alpha() * N.y();
  // div-div
  //  dst(1,0) = ST.alpha() * N.x();
  //  dst(2,0) = ST.alpha() * N.y();
}

/*-----------------------------------------*/

void StokesGls2d::LMatrix(nmatrix<double>& A, 
			const FemFunction& U,
			const TestFunction& V) const
{
  A(0,1) = V.x();
  A(0,2) = V.y();
  A(1,0) = V.x();
  A(2,0) = V.y();
}
}

/*-----------------------------------------*/

