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


#ifndef __gmresclass_h
#define __gmresclass_h

#include "mult.h"
#include "cginfo.h"
#include "nvector.h"

/********************************************************************/

namespace Gascoigne
{
template<class SOLVER, class PRECONDITIONER, class VECTOR>
class GMRES
{
  typedef nvector<double>  dvector;

  nvector<dvector> H;
  nvector<double>  gamma, ci, si;

  std::vector<VECTOR> mem;
  
  int vmax, left_precondition;
  
  void   new_memory          ();
  void   givens_rotation     (dvector&, int);
  void   solution            (VECTOR&, VECTOR&, int);
  double orthogonalization   (dvector&, int, VECTOR&) const;
  bool   reortho_test        (const VECTOR&, double) const;

  SOLVER&         A;
  PRECONDITIONER& P;
  
 public:

  GMRES(SOLVER&, PRECONDITIONER&, int);
  ~GMRES();
  void init();

  int solve          (VECTOR& x, const VECTOR& b, CGInfo& info);
  int restarted_solve(VECTOR& x, const VECTOR& b, CGInfo& info);
};
}

/********************************************************************/

#endif
