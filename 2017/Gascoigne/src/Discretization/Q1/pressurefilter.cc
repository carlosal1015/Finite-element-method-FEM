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


#include "pressurefilter.h"

/*-----------------------------------------*/

namespace Gascoigne
{
PressureFilter::PressureFilter() : nvector<double>(), domainsize(0.),  nhanging(0) {}

/*-----------------------------------------*/

PressureFilter::~PressureFilter() {}

/*-----------------------------------------*/

void PressureFilter::ReInit(int n, int nhn)
{
  resize(n);
  zero();
  domainsize = 0.;
  nhanging = nhn;
}

/*-----------------------------------------*/

DoubleVector PressureFilter::IntegrateVector(const GlobalVector& u) const
{
  assert(size());
  assert(Active());
  
  DoubleVector dst(u.ncomp(),0.);
  
  for (int i=0; i<component.size(); i++)
    {
      int c = component[i];
      for (int j=0; j<u.n(); j++)
	{
	  dst[c] += u(j,c)* (*this)[j];
	}      
    }
  return dst;
}

/*-----------------------------------------*/

void PressureFilter::SubtractMean(GlobalVector& u) const
{
  assert(size());
  assert(domainsize>0.);
  DoubleVector mean = IntegrateVector(u);
  
  for (int i=0; i<component.size(); i++)
    {
      int   comp = component[i];
      double sub = mean[comp]/domainsize;
      u.CompAdd(comp,-sub);
    }  
}

/*-----------------------------------------*/

void PressureFilter::SubtractMeanAlgebraic(GlobalVector& u) const
{
  for (int i=0; i<component.size(); i++)
    {
      int comp = component[i];
      double d = 0.;
      for (int j=0; j<u.n(); j++)
	{
	  d += u(j,comp);
	}      
      d /= u.n() - nhanging;
      u.CompAdd(comp,-d);
    }
}

/*-----------------------------------------*/

}

