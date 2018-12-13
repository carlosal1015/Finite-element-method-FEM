/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#include  "navierstokeslps2d.h"
#include  "filescanner.h"

namespace Gascoigne
{

/*-----------------------------------------*/

NavierStokesLps2d::NavierStokesLps2d(const ParamFile* filename) 
  : LpsEquation(), NavierStokes2d() 
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("cut"  , &_cut  ,  1.e8);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("delta", &ST.delta0, 0.25);
  DFH.insert("tau",   &ST.tau0,  0.);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("penalty",&_penalty, 0.);
  DFH.insert("dtfactor",&ST.dtfactor(), 1.); // denke das klappt nur mit dtfactor=0 ... aber 1 ist bei ST default.. :/

  FileScanner FS(DFH,filename,"Equation");
}

/*-----------------------------------------*/

void NavierStokesLps2d::lpspoint(double h, const FemFunction& U, const Vertex2d& v)const
{
  ST.ReInit(h,_visc,U[1].m(),U[2].m());
}
 
/*-----------------------------------------*/

void NavierStokesLps2d::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
{
  b[0] += ST.alpha() * (UP[0].x()*N.x()+UP[0].y()*N.y());
  double betaN = Convection(U,N);
  double betaU1 = U[1].m()*UP[1].x()+U[2].m()*UP[1].y();
  double betaU2 = U[1].m()*UP[2].x()+U[2].m()*UP[2].y();
  b[1] += ST.delta() * betaU1*betaN;
  b[2] += ST.delta() * betaU2*betaN;
  //b[1] += ST.delta() * (UP[1].x()*N.x()+UP[1].y()*N.y());
  //b[2] += ST.delta() * (UP[2].x()*N.x()+UP[2].y()*N.y());

  double divU = Divergence(UP);
  b[1] += ST.tau() * divU * N.x();
  b[2] += ST.tau() * divU * N.y();
}
 
/*-----------------------------------------*/

void NavierStokesLps2d::StabMatrix(EntryMatrix& A,  const FemFunction& U, 
 const TestFunction& Np, const TestFunction& Mp) const
{
  double laplace = Laplace(Mp,Np);
  double betaM = Convection(U,Mp);
  double betaN = Convection(U,Np);
  double betabeta = ST.delta() * betaM*betaN;
  //double betabeta = ST.delta() * laplace;

  A(0,0) += ST.alpha() * laplace;
  A(1,1) += betabeta;
  A(2,2) += betabeta;

  A(1,1) += ST.tau() * Mp.x() * Np.x();
  A(1,2) += ST.tau() * Mp.y() * Np.x();
  A(2,1) += ST.tau() * Mp.x() * Np.y();
  A(2,2) += ST.tau() * Mp.y() * Np.y();
}
}
