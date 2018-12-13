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


#ifndef __Q43d_h
#define __Q43d_h

#include "q4.h"

namespace Gascoigne
{

/**********************************************************/

  class Q43d : public Q4
  {
    protected:
      int GetPatchNumber(const Vertex3d& p0, Vertex3d& p) const;
      void VertexTransformation(const Vertex3d& p0, Vertex3d& p, int iq) const;

    public:
      Q43d();
      ~Q43d();

      std::string GetName() const {return "Q43d";}

      void BasicInit(const ParamFile* paramfile);
  };

/**********************************************************/

}

#endif
