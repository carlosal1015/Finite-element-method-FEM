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


#ifndef  __waveequation_h
#define  __waveequation_h

#include  "equation.h"
#include  "boundaryequation.h"

/*---------------------------------------------------*/

class WaveEquation : public Gascoigne::Equation
{
  double c2;
  double Laplace(const Gascoigne::TestFunction& U, const Gascoigne::TestFunction& N) const;

public:

  WaveEquation();

  int  GetNcomp() const { return 1; }

  std::string GetName() const { return "WaveEquation";}

  void SetTimePattern(Gascoigne::TimePattern& TP) const 
  { TP(0,0) = 1.;}

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
	    const Gascoigne::TestFunction& N) const;

  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
	      const Gascoigne::TestFunction& M, 
	      const Gascoigne::TestFunction& N) const;
};

/*---------------------------------------------------*/

class WaveBoundaryEquation : public Gascoigne::BoundaryEquation
{
  // for absorbing boundary condition
  //
  double c2;

public:

  WaveBoundaryEquation();
  int GetNcomp()const { return 1;}

  std::string GetName() const { return "WaveBoundaryEquation";}

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
	    const Gascoigne::TestFunction& N, int color) const;

  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
	      const Gascoigne::TestFunction& M,
              const Gascoigne::TestFunction& N, int color) const;
};

/*---------------------------------------------------*/

#endif
