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


#include "edgeinfo.h"

using namespace std;

/**********************************************************/

namespace Gascoigne
{
template<int DIM>
void EdgeInfo<DIM>::BasicInit(const Edge* edge, int ncomp, const fixarray<2*DIM-2,int>& vertex)
{
  _count  = 0;
  _edge   = edge;
  _vertex = vertex;
  _u.ReInit(ncomp,2*DIM-2);
  _u.zero();
}

/**********************************************************/

template<int DIM>
void EdgeInfo<DIM>::AddNodes(const LocalVector& u)
{
  for (int i=0; i<2*DIM-2; i++)
    {
      for (int c=0; c<_u.ncomp(); c++)
	{
	  _u(i,c) += u(i,c);
	}
    }
  _count++;
}

/**********************************************************/

template<int DIM>
fixarray<2*DIM-2,double> EdgeInfo<DIM>::GetNorm() const
{
  fixarray<2*DIM-2,double> norm(0.);

  for (int i=0; i<2*DIM-2; i++)
    {
      for (int c=0; c<_u.ncomp(); c++)
	{
	  norm[i] += _u(i,c) * _u(i,c);
	}
    }
  return norm;
}

/**********************************************************/

template<int DIM>
void EdgeInfo<DIM>::ShowStatistics() const
{
  cout << _count << ": " << _vertex << " = " << _u << endl;
}

/**********************************************************/

template class EdgeInfo<2>;
template class EdgeInfo<3>;
}
