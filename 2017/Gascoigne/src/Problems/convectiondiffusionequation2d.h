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


#ifndef __ConvectionDiffusionEquation_h
#define __ConvectionDiffusionEquation_h

#include "equation.h"
#include "paramfile.h"

/*---------------------------------------------------*/

namespace Gascoigne{

class ConvectionDiffusionEquation2d : public virtual Equation
{
protected:

  double _visc, _bx, _by;

  double betax() const {return _bx;}
  double betay() const {return _by;}
  double Convection(const TestFunction& N) const;

public:

  ConvectionDiffusionEquation2d(const ParamFile* paramfile);
  ~ConvectionDiffusionEquation2d() {}
  void SetTimePattern(TimePattern& P) const;

  std::string GetName() const {return "ConvectionDiffusionEquation";}

  int  GetNcomp()    const { return 1;}

  void OperatorStrong(DoubleVector& b, const FemFunction& U) const;

  void point(double h, const FemFunction& U, FemData& Q, 
	     const Vertex2d& v) const {}

  void Form(VectorIterator b, const FemFunction& U, 
	    const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, 
	      const TestFunction& M, const TestFunction& N) const;
};
}
/*---------------------------------------------------*/

#endif
