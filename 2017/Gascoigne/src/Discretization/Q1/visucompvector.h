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


#ifndef __visucompvector_h
#define __visucompvector_h

#include  "compvector.h"
#include  "visudata.h"

/***************************************************************/

namespace Gascoigne
{
class VisuCompVector : public VisuData
{
 protected:

  typedef DoubleVector::iterator         pointer;
  typedef DoubleVector::const_iterator   const_pointer;
  
  const DoubleVector *uR, *zR;
  int           ncomp;

 public:
  
  VisuCompVector(const DoubleVector& u) : uR(&u), zR(0)
    {
      ncomp = u.ncomp();
    }
  VisuCompVector(const DoubleVector& u, const DoubleVector& z) : uR(&u), zR(&z) 
    {
      ncomp = u.ncomp()+z.ncomp();
    }

  int    visucomp()            const { return ncomp;}
  int    visun()               const { return uR->n();}
  double visudata(int i,int c) const 
    { 
      if ( c < uR->ncomp() ) return (*uR)(i,c);
      return (*zR)(i,c-uR->ncomp());
    }
};
}

/***************************************************************/

#endif
