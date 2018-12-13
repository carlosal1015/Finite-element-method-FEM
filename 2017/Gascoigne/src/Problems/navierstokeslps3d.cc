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


#include  "navierstokeslps3d.h"
#include  "filescanner.h"

using namespace std;

namespace Gascoigne
{

/*-----------------------------------------*/

NavierStokesLps3d::NavierStokesLps3d(const ParamFile* filename) 
  : NavierStokes3d()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("cut"  , &_cut  ,  1.e8);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("delta", &ST.delta0, 0.25);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH,filename,"Equation");
}

/*-----------------------------------------*/

string NavierStokesLps3d::GetName() const
{
  return "NavierStokesLps3d";
}

/*-----------------------------------------*/

void NavierStokesLps3d::lpspoint(double h, const FemFunction& U, const Vertex3d& v)const
{
  ST.ReInit(h,_visc,U[1].m(),U[2].m(),U[3].m());
}
 
/*-----------------------------------------*/

void NavierStokesLps3d::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
{
  b[0] += ST.alpha() * Laplace(UP[0],N);
  double betaN  = Convection(U,N);
  double betaU1 = Convection(U,UP[1]);
  double betaU2 = Convection(U,UP[2]);
  double betaU3 = Convection(U,UP[3]);
  b[1] += ST.delta() * betaU1*betaN;
  b[2] += ST.delta() * betaU2*betaN;
  b[3] += ST.delta() * betaU3*betaN;
}
 
/*-----------------------------------------*/

void NavierStokesLps3d::StabMatrix(EntryMatrix& A,  const FemFunction& U, 
 const TestFunction& Np, const TestFunction& Mp) const
{
  double laplace = Laplace(Mp,Np);
  double betaM = Convection(U,Mp);
  double betaN = Convection(U,Np);
  double betabeta = ST.delta() * betaM*betaN;

  ////////////// Continuity ////////////////////////////////////////////////

  A(0,0) += ST.alpha() * laplace;
  A(1,1) += betabeta;
  A(2,2) += betabeta;
  A(3,3) += betabeta;
}
}

/*-----------------------------------------*/
