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


#ifndef  __NavierStokesLps3d_h
#define  __NavierStokesLps3d_h

#include  "navierstokes3d.h"
#include  "lpsequation.h"
#include  "lpsstabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class NavierStokesLps3d : public NavierStokes3d, public virtual LpsEquation
{
  protected:

  mutable LpsStabilization ST;
  
  public:

  NavierStokesLps3d(const Gascoigne::ParamFile* filename);

  std::string GetName() const;

  void SetTime(double time, double dt) const {Application::SetTime(time,dt); ST.DeltaT() = dt;}
  //
  /// Computation of lps stabilization parameters
  //
  void lpspoint(double h, const FemFunction& U, const Vertex3d& v) const;
  //
  /// for local-projection stabilization (lps)
  //
  void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
    
  void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;
};
}

#endif
