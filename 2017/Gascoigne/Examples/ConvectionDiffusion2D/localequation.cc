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


#include  "localequation.h"
#include  "filescanner.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

LocalEquation::LocalEquation(const ParamFile* paramfile) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("sigma",&sigma,0.);
  DFH.insert("visc",&visc,1.);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  FileScanner FS(DFH,paramfile,"Equation");
}

/* ----------------------------------------- */

void LocalEquation::glspoint(double h, const FemFunction& U, const Vertex2d& v) const
{
  ST.ReInit(h,visc,betax(),betay());
}

/* ----------------------------------------- */

void LocalEquation::SetFemData(FemData& Q) const
{
  assert(Q.count("beta"));
  q = &Q["beta"];
}

/* ----------------------------------------- */

void LocalEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  b[0] += visc* (U[0].x()*N.x()+U[0].y()*N.y());
  b[0] += betax()* U[0].x()*N.m();
  b[0] += betay()* U[0].y()*N.m();
}

/* ----------------------------------------- */

void LocalEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  A(0,0) += visc* (M.x()*N.x()+M.y()*N.y());
  A(0,0) += betax()* M.x()*N.m();
  A(0,0) += betay()* M.y()*N.m();
}

/* ----------------------------------------- */

void LocalEquation::L(nvector<double>& dst, const FemFunction& U) const
{
  dst[0] = betax()* U[0].x() + betay()* U[0].y();
}

/* ----------------------------------------- */

void LocalEquation::S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const
{
  dst(0,0) = ST.alpha()*(betax()* N.x() + betay()* N.y());
}

/* ----------------------------------------- */

void LocalEquation::LMatrix(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const
{
  dst(0,0) = betax()* N.x() + betay()* N.y();
}
