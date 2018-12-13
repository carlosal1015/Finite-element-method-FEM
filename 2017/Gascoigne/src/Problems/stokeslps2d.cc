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


#include  "stokeslps2d.h"
#include  "filescanner.h"

/*-----------------------------------------*/

namespace Gascoigne
{
StokesLps2d::~StokesLps2d()
{
}

/*-----------------------------------------*/

StokesLps2d::StokesLps2d() :  LpsEquation(), Stokes2d()
{
  _penalty = 0.; _visc = 1.;
  
  ST.xeta0 = 6.;
  ST.alpha0 = 0.2;
}
 
/*-----------------------------------------*/

StokesLps2d::StokesLps2d(const ParamFile* filename) : 
  LpsEquation(), Stokes2d()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("xdtfactor" , &ST.dtfactor(), 1.);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH,filename,"Equation");
}

/*-----------------------------------------*/

void StokesLps2d::lpspoint(double h, const FemFunction& U, const Vertex2d& v)const
{
//   ST.alpha() = 0.1 * (h*h)/_visc;
  ST.ReInit(h,_visc);
}

/*-----------------------------------------*/

void StokesLps2d::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
{
  b[0] += ST.alpha() * (UP[0].x()*N.x()+UP[0].y()*N.y());
}
 
/*-----------------------------------------*/

void StokesLps2d::StabMatrix(EntryMatrix& A,  const FemFunction& U, 
 const TestFunction& Np, const TestFunction& Mp) const
{
  double laplace = Laplace(Mp,Np);

  A(0,0) += ST.alpha() * laplace;
}
}
