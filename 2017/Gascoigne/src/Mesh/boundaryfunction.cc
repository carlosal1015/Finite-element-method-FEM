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


#include  "boundaryfunction.h"

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne
{
template<int DIM>
void BoundaryFunction<DIM>::grad(Vector& dst, const Vector& src) const 
{
  double eps = 1e-6;
  
  for(int i=0;i<DIM;i++)
    {
      Vector cl(src), cr(src);
      cl[i] -= eps;
      cr[i] += eps;
      dst[i] = ((*this)(cr)-(*this)(cl))/(2.*eps);
    }
}

/*---------------------------------------------------*/

template<int DIM>
void BoundaryFunction<DIM>::newton(Vector& dst) const
{
  int    maxi = 10;
  double tol  = 1.e-12;

  Vector z;

  grad(z,dst);
  double res = (*this)(dst);

  for (int i=0; (i<maxi) && (fabs(res)>tol) ; i++)
    {
      Vector zz;
      grad(zz,dst);
      double bgrad = z*zz;

      if (fabs(bgrad)<=1.e-15)
	{
	  cerr << "BoundaryFunction<DIM>::newton()\n";
	  cerr << "Grad=0 in boundary_newton (res= "<< res << " )\n";
	  cerr << "iter " << i << endl;
	  abort();
	}
      res /= bgrad;
      dst.add(-res,z);
      res = (*this)(dst);
    }
  if (fabs(res)>tol)
    {
      cerr << "BoundaryFunction<DIM>::newton()\n";
      cerr << "No Convergence in boundary_newton (res= "<< res << " )\n";
      abort();
    }
}

/*---------------------------------------------------*/

template class BoundaryFunction<2>;
template class BoundaryFunction<3>;
}
