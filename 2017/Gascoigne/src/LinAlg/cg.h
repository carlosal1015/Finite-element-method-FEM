/**
*
* Copyright (C) 2004, 2008 by the Gascoigne 3D authors
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


#ifndef __cg_h 
#define __cg_h 

#include "cginfo.h" 

/********************************************************************/

// S = Solver
// V = Vector

namespace Gascoigne
{
template <class S, class V>
class CG
{
  S   &solver;

public:

  CG(S& s) : solver(s) {}

  void solve(V& x, const V& b, CGInfo& info);
};
}

/********************************************************************/


#endif 
