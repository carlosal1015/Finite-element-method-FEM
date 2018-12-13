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


#ifndef __numderivative_h
#define __numderivative_h

namespace Gascoigne
{
template <class C, class MAT, class VEC>
void numderivative(C& application, MAT& M, const VEC& x, double eps=1.e-4)
{
  int m = x.size();

  VEC up(x), u(x), xp(x);
  up.zero(); u.zero();

//  double ieps = 1./eps;
  
  application.f(u,x);
  
  for (int i=0; i<m ; i++)
    {
      //xp[i] += eps;
      xp[i] *= 1.+eps;
      application.f(up,xp);
      for (int j=0; j<m ; j++)
	{
	  //M(j,i) = (up[j]-u[j])*ieps;
	  M(j,i) = (up[j]-u[j])/(xp[i]*eps);
	}
      //xp[i] -= eps;
      xp[i] /= 1.+eps;
    }
}
}

#endif
