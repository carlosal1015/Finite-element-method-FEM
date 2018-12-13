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


#ifndef  __StokesLps3d_h
#define  __StokesLps3d_h

#include  "stokes3d.h"
#include  "lpsequation.h"
#include  "stabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class StokesLps3d : public Stokes3d, public virtual LpsEquation
{
protected:

  //
  /// handles the stabilization parameters
  //
  mutable Stabilization ST;

public:

  ~StokesLps3d();
  StokesLps3d();
  StokesLps3d(const ParamFile* filename);

  std::string GetName() const { return "StokesLps3d";}

  void SetTime(double time, double dt) const {Application::SetTime(time,dt); ST.DeltaT() = dt;}
  //
  /// Computation of lps stabilization parameters
  //
  void lpspoint(double h, const FemFunction& U, const Vertex3d& v) const;

  //
  /// Stablization terms
  //
  void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
    
  void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;
};
}

#endif
