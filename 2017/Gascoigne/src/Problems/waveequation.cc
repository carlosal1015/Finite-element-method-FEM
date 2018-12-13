/**
*
* Copyright (C) 2010 by the Gascoigne 3D authors
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


#include  "waveequation.h"
#include  "filescanner.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

WaveEquation::WaveEquation() : 
  Equation() { c2=1.;}

/*---------------------------------------------------*/

double WaveEquation::Laplace(const TestFunction& U,
			     const TestFunction& N) const
{
  return U.x()*N.x() + U.y()*N.y();
}

/*---------------------------------------------------*/

void WaveEquation::Form(VectorIterator b, const FemFunction& U, 
			const TestFunction& N) const
{
  b[0] += c2 * Laplace(U[0],N);
}

/*---------------------------------------------------*/

void WaveEquation::Matrix(EntryMatrix& A, const FemFunction& U, 
			  const TestFunction& M, const TestFunction& N) const
{
  A(0,0) += c2 * Laplace(M,N);
}

/*-----------------------------------------*/

WaveBoundaryEquation::WaveBoundaryEquation() : 
  BoundaryEquation()  { c2=1.;}

/*---------------------------------------------------*/

void WaveBoundaryEquation::Form(VectorIterator b, const FemFunction& U, 
				const TestFunction& N, int color) const
{ 
  b[0] += c2*U[0].m()*N.m();
}

/*---------------------------------------------------*/

void WaveBoundaryEquation::Matrix(Gascoigne::EntryMatrix& A, 
				  const Gascoigne::FemFunction& U, 
				  const Gascoigne::TestFunction& M,
				  const Gascoigne::TestFunction& N, int color) const
{
  A(0,0) += c2*M.m()*N.m();
}

/*---------------------------------------------------*/
