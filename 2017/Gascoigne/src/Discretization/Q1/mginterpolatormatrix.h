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


#ifndef  __MgInterpolatorMatrix_h
#define  __MgInterpolatorMatrix_h

#include  "mginterpolatorinterface.h"
#include  "columnstencil.h"
#include  "gascoigne.h"


/*-----------------------------------------*/


namespace Gascoigne
{
class MgInterpolatorMatrix : public virtual MgInterpolatorInterface
{
private:

  ColumnStencil  ST;
  DoubleVector      val;

public:


  MgInterpolatorMatrix() : MgInterpolatorInterface() {}

  ColumnStencil& GetStencil() {return  ST;}
  const ColumnStencil& GetStencil() const {return  ST;}

  DoubleVector& GetAlpha() {return val;}
  double Alpha(int pos) const {return val[pos];}

  void restrict_zero   (GlobalVector& uL, const GlobalVector& ul) const;
  void prolongate_add  (GlobalVector& ul, const GlobalVector& uL) const;
  void SolutionTransfer(GlobalVector& uL, const GlobalVector& ul) const;

};
}

#endif
