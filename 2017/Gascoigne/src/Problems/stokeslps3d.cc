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


#include "stokeslps3d.h"
#include "filescanner.h"

/*-----------------------------------------*/

namespace Gascoigne
{
StokesLps3d::~StokesLps3d() {}

/*-----------------------------------------*/

StokesLps3d::StokesLps3d() : 
  LpsEquation(), Stokes3d() 
{ 
  _penalty = 0.; _visc = 1.;
  
  ST.xeta0 = 6.;
  ST.alpha0 = 0.2;
}

/*-----------------------------------------*/

StokesLps3d::StokesLps3d(const ParamFile* filename) : 
  LpsEquation(), Stokes3d() 
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH,filename,"Equation");
}

/*-----------------------------------------*/

void StokesLps3d::lpspoint(double h, const FemFunction& U, const Vertex3d& v)const
{
  ST.ReInit(h,_visc);
}

/*-----------------------------------------*/

void StokesLps3d::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
{
  b[0] += ST.alpha() * Laplace(UP[0],N);
}
 
/*-----------------------------------------*/

void StokesLps3d::StabMatrix(EntryMatrix& A,  const FemFunction& U, const TestFunction& Np, 
			     const TestFunction& Mp) const
{
  A(0,0) += ST.alpha() * Laplace(Mp,Np);
}
}
