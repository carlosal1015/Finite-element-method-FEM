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


#ifndef __gascoignemath_h
#define __gascoignemath_h

#include "math.h"

namespace Gascoigne
{
inline double pi()
{
  return 3.14159265358979323846;
}

inline double max(double a, double b) 
{
  if (a>b) return a;
  return b;
}

inline int max_int(int a, int b) 
{
  if (a>b) return a;
  return b;
}

inline double min(double a, double b) 
{
  if (a<b) return a;
  return b;
}

inline int min_int(int a, int b) 
{
  if (a<b) return a;
  return b;
}

inline int abs_int(int a)
{
  if (a>0) return a;
  return -a;
}

/* #define PI 3.14159265358979323846 */

/* #ifndef MIN */
/* #define MIN(a,b) ( ((a)>(b)) ? (b) : (a) ) */
/* #endif */

/* #ifndef ABS */
/* #define ABS(a)   ( ((a)>(0)) ? (a) : (-a) ) */
/* #endif */

/* #ifndef FRAC */
/* #define FRAC(x)  (fabs(x-static_cast<int>(x))) */
/* #endif */

/* int ggt(int n1, int n2); */
}

#endif

