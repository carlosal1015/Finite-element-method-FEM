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


#ifndef  __NavierStokesSplit_h
#define  __NavierStokesSplit_h

#include  "lpsequation.h"
#include  "lpsstabilization.h"
#include  "paramfile.h"

/*---------------------------------------------------*/

class NavierStokesSplitLps2d : public virtual Gascoigne::LpsEquation
{
protected:

  mutable double _h, _visc;
  mutable Gascoigne::LpsStabilization ST;

  double Laplace(const Gascoigne::TestFunction& U, const Gascoigne::TestFunction& N) const;
  double Convection(const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;
  double Divergence(const Gascoigne::FemFunction& U) const;

public:

  NavierStokesSplitLps2d(const Gascoigne::ParamFile* filename);

  std::string GetName() const { return "NavierStokesSplitLps2d";}
  int  GetNcomp()       const { return 2; }

  void SetTime(double time, double dt) const;

  void SetTimePattern(Gascoigne::TimePattern& P) const;

  void point(double h, const Gascoigne::FemFunction& U, const Gascoigne::Vertex2d& v) const 
  { _h = h;}

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
	    const Gascoigne::TestFunction& N) const;

  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
	      const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;
  //
  /// Computation of lps stabilization parameters
     //
    void lpspoint(double h, const Gascoigne::FemFunction& U, 
		  const Gascoigne::Vertex2d& v) const
  {
    ST.ReInit(h,_visc,U[0].m(),U[1].m());
  }
  //
  /// for local-projection stabilization (lps)
  //

  void StabForm(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
		const Gascoigne::FemFunction& UP, const Gascoigne::TestFunction& N) const;
    
  void StabMatrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
		  const Gascoigne::TestFunction& Np, const Gascoigne::TestFunction& Mp) const;
};

/*---------------------------------------------------*/

#endif
