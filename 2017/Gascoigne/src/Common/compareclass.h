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


#ifndef  __compareclass_h
#define  __compareclass_h

#include "nvector.h"
#include <algorithm>

/*********************************************************/

namespace Gascoigne
{
template<class T>
class CompareLess
{
  protected:
  
  const T*   wp;

  public:

  CompareLess<T>() {}
  CompareLess<T>(const T& w)
    { 
      wp  = &w;
    }
  bool operator()(int i1,int i2) const
    {
      if( (*wp)[i1]  < (*wp)[i2] ) return 1;
      return 0;
    }
};

/*********************************************************/

template<class T>
class CompareObject
{
  protected:
  
  const T& P;

  public:

  CompareObject(const T& p) : P(p)
    {};
  
  bool operator()(int i, int j) const
    { return P[i]<P[j]; }
};

/*********************************************************/

template<class T>
class CompareObjectBigToSmall
{
  protected:
  
  const T& P;
  
  public:
  
  CompareObjectBigToSmall(const T& p) : P(p)
    {};
  
  bool operator()(int i, int j) const
    { return P[j]<P[i]; }
};
}

#endif

