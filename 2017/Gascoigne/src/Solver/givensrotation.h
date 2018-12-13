/**
*
* Copyright (C) 2010 by the Gascoigne 3D authors
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


#ifndef __givensrotation_h
#define __givensrotation_h

#include "solverinterface.h"

/*---------------------------------------------------------------*/
 
namespace Gascoigne
{
class GivensRotation
{
  //
  // Data
  //
  int                  n;
  DoubleMatrix  H;
  DoubleVector  ci, si,gamma;

 public:
  
  GivensRotation(int nn, double) ;
  double&         matrix(int i, int j)             { return H(i,j);}
  void               givens(double& a, double& b, int i) const;
  double           orthogonalization(int dim) ;
  DoubleVector getcoefficients();
};
} 

/*---------------------------------------------------------------*/
 
#endif
