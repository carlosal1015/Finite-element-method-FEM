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


#ifndef  __PressureFilter_h
#define  __PressureFilter_h

#include "nvector.h"
#include "gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class PressureFilter : public DoubleVector 
{
 protected:

  IntVector component;
  double       domainsize;
  int          nhanging;

 public:

  PressureFilter() ;
  ~PressureFilter();

  void SetComponents(const IntVector& c) { component = c;}
  bool Active() const { return component.size()>0;}

  void ReInit(int n, int nhn);

  void AddDomainPiece(double val) { domainsize += val;}

  DoubleVector IntegrateVector(const GlobalVector& u) const;
  void SubtractMean(GlobalVector& u) const;
  void SubtractMeanAlgebraic(GlobalVector& u) const;
};
}

/*-----------------------------------------*/

#endif
