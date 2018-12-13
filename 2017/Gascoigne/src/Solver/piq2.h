/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#ifndef __PiQ2_h
#define __PiQ2_h

#include "gascoignemesh.h"

namespace Gascoigne
{

/**********************************************************/

  class PiQ2
  {
    protected:
      const GascoigneMesh   *_MP;
      nvector<DoubleVector>  _q2weight;

    public:
      PiQ2();
      ~PiQ2() { }

      void Init(const MeshInterface* MI);
      void vmult(GlobalVector& y, const GlobalVector& x) const;
  };

/**********************************************************/

}

#endif
