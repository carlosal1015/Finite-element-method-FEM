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


#ifndef __visudatascalar_h
#define __visudatascalar_h

#include  "visudata.h"
#include  "nvector.h"

/***************************************************************/

namespace Gascoigne
{
class VisuDataScalar : public VisuData
{
 protected:

  const DoubleVector& vR;

 public:

  virtual ~VisuDataScalar(){}
  VisuDataScalar(const DoubleVector& v) : vR(v) {}

  virtual int    visucomp()     const {return 1;}
  virtual int    visun()        const {return vR.size();}
  virtual double visudata(int i,int c) const { return vR[i];}
};
}

/***************************************************************/

#endif
