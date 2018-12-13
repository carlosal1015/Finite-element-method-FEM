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


#ifndef __NewmarkSolver_h
#define __NewmarkSolver_h

#include  "timesolver.h"

/*----------------------------------------------------------------------------*/

class NewmarkSolver : public Gascoigne::TimeSolver
{
public:

  NewmarkSolver() : Gascoigne::TimeSolver() {}
  void Form(Gascoigne::VectorInterface& gy, const Gascoigne::VectorInterface& gx, double d) const;
  void FormWithoutMass(Gascoigne::VectorInterface& gy, const Gascoigne::VectorInterface& gx, 
		       double d, double s=1.) const;
  void AssembleMatrix(const Gascoigne::VectorInterface& gu, double d);
};

/*----------------------------------------------------------------------------*/

#endif
