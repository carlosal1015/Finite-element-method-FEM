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


#ifndef __mult_h
#define __mult_h

#include "gascoignemath.h"

namespace Gascoigne
{
template<class MATRIX1, class MATRIX2, class MATRIX3>
void multeq(MATRIX1& M, const MATRIX2& A, const MATRIX3& B)
{
  typename MATRIX1::iterator         p  = M.begin();
  typename MATRIX2::const_iterator   Ap = A.begin();
  typename MATRIX3::const_iterator   Bp;
  
  int m  = M.m();
  int ma = A.m();
  int mb = B.m();
  int nb = B.n();
  
  while(p!=M.end())
    {
      for(int j=0;j<m;j++)
	{
	  Bp = B.begin()+j;
	  *p = 0.;
	  for(int k=0;k<nb;k++)
	    {
	      *p += (*Ap++) * (*Bp);
	      Bp += mb;
	    }
	  p++;
	  Ap -= nb;
	}
      Ap += ma;
    }
}

template<class MATRIX, class VECTOR>
void vmult(VECTOR& y, const MATRIX& M, const VECTOR& x)
{
  typename VECTOR::iterator         py = y.begin();
  typename VECTOR::const_iterator   px;
  typename MATRIX::const_iterator   p = M.begin();

  while(p!=M.end())
    {
      px = x.begin();
      for(int j=0;j<M.m();j++)
	{
	  *py += (*p++) * (*px++);
	}
      py++;
    }
}

template<class MATRIX, class VECTOR1, class VECTOR2>
void vmulteq(VECTOR1& y, const MATRIX& M, const VECTOR2& x)
{
  typename VECTOR1::iterator         py = y.begin();
  typename VECTOR2::const_iterator   px;
  typename MATRIX::const_iterator   p = M.begin();

  while(p!=M.end())
    {
      px = x.begin();
      *py = 0.;
      for(int j=0;j<M.m();j++)
	{
	  *py += (*p++) * (*px++);
	}
      py++;
    }
}

template<class MATRIX, class ITERATOR, class CITERATOR>
void vmulteq2(ITERATOR& py, const MATRIX& M, CITERATOR& px0)
{
  CITERATOR                         px;
  typename MATRIX::const_iterator   p = M.begin();

  while(p!=M.end())
    {
      px = px0;
      *py = 0.;
      for(int j=0;j<M.m();j++)
	{
	  *py += (*p++) * (*px++);
	}
      py++;
    }
}

template<class MATRIX, class CITERATOR>
double norme(const MATRIX& M, CITERATOR px0)
{
  double d = 0.;
  typename MATRIX::const_iterator   p = M.begin();

  CITERATOR px,py=px0;
  while(p!=M.end())
    {
      px = px0;
      for(int j=0;j<M.m();j++)
	{
	  d += (*py) * (*p++) * (*px++);
	}
      py++;
    }
  return d;
}

template<class MATRIX, class VECTOR>
void backward(VECTOR& dst, const MATRIX& M, const VECTOR& src)
    {
      int i,j;
      int nu = min_int(M.m(),M.n());
      double s;
      for (i=nu-1;i>=0;i--)
	{
	  s = src[i];
	  for (j=i+1;j<nu;j++) s -= dst[j] * M(i,j);
	  dst[i] = s/M(i,i);
	}
    }
}

#endif

