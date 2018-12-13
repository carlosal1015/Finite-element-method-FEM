/**
*
* Copyright (C) 2008 by the Gascoigne 3D authors
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


#include  "navierstokessplit.h"
#include  "lpsstabilization.h"
#include  "filescanner.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

double NavierStokesSplitLps2d::Laplace(const TestFunction& U, 
				       const TestFunction& N) const
{
  return U.x()*N.x() + U.y()*N.y();
}

/*---------------------------------------------------*/

double NavierStokesSplitLps2d::Convection(const FemFunction& U, 
					  const TestFunction& N) const
{
  return U[0].m()*N.x() + U[1].m()*N.y();
}

/*---------------------------------------------------*/

double NavierStokesSplitLps2d::Divergence(const FemFunction& U) const
{
  return U[0].x() + U[1].y();
}

/*---------------------------------------------------*/

NavierStokesSplitLps2d::NavierStokesSplitLps2d(const ParamFile* pf) : 
  LpsEquation()
{
  _h = 0.;
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 1.);
  DFH.insert("delta", &ST.delta0, 0.25);
  DFH.insert("tau",   &ST.tau0,  0.);
  DFH.insert("xeta" , &ST.xeta0, 6.);
      
  FileScanner FS(DFH, pf, "Equation");
}

/*---------------------------------------------------*/

void NavierStokesSplitLps2d::SetTime(double time, double dt) const 
{ 
  Application::SetTime(time,dt); 
  ST.DeltaT() = dt;
}

/*---------------------------------------------------*/

void NavierStokesSplitLps2d::SetTimePattern(TimePattern& P) const
{
  P(0,0) = 1.;
  P(1,1) = 1.;
}

/*---------------------------------------------------*/

void NavierStokesSplitLps2d::Form(VectorIterator b, const FemFunction& U, 
				  const TestFunction& N) const
{
  ////////////// Momentum ////////////////////////////////////////////////
  
  b[0] += Convection(U,U[0]) * N.m();
  b[1] += Convection(U,U[1]) * N.m();
  
  // viscous terms
  b[0] += _visc * Laplace(U[0],N);
  b[1] += _visc * Laplace(U[1],N);
}

/*---------------------------------------------------*/

void NavierStokesSplitLps2d::Matrix(EntryMatrix& A, 
				    const FemFunction& U, 
				    const TestFunction& M, 
				    const TestFunction& N) const
{
  double MN = M.m()*N.m();
  double laplace = Laplace(M,N);
  
  ////////////// Momentum ////////////////////////////////////////////////
  
  double cl = Convection(U,M) * N.m() + _visc*laplace;
  
  A(0,0) += cl;
  A(1,1) += cl;
  
  A(0,0) += U[0].x()*MN;
  A(0,1) += U[0].y()*MN;
  A(1,0) += U[1].x()*MN;
  A(1,1) += U[1].y()*MN;
}

/*---------------------------------------------------*/

void NavierStokesSplitLps2d::StabForm(VectorIterator b, const FemFunction& U, 
				      const FemFunction& UP, const TestFunction& N) const
{
  double betaN = Convection(U,N);
  double betaU1 = Convection(U,UP[0]);
  double betaU2 = Convection(U,UP[1]);
  b[0] += ST.delta() * betaU1*betaN;
  b[1] += ST.delta() * betaU2*betaN;
  
  double divU = Divergence(UP);
  b[0] += ST.tau() * divU * N.x();
  b[1] += ST.tau() * divU * N.y();
}

/*---------------------------------------------------*/
    
void NavierStokesSplitLps2d::StabMatrix(EntryMatrix& A, const FemFunction& U, 
					const TestFunction& Np, const TestFunction& Mp) const
{
  double betaM = Convection(U,Mp);
  double betaN = Convection(U,Np);
  double betabeta = ST.delta() * betaM*betaN;
  
  A(0,0) += betabeta;
  A(1,1) += betabeta;
  
  A(0,0) += ST.tau() * Mp.x() * Np.x();
  A(0,1) += ST.tau() * Mp.y() * Np.x();
  A(1,0) += ST.tau() * Mp.x() * Np.y();
  A(1,1) += ST.tau() * Mp.y() * Np.y();
}

/*---------------------------------------------------*/
