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


#ifndef __Q42d_h
#define __Q42d_h

#include "q4.h"

namespace Gascoigne
{

/**********************************************************/

  class Q42d : public Q4
  {
    protected:
      int GetPatchNumber(const Vertex2d& p0, Vertex2d& p) const;
      void VertexTransformation(const Vertex2d& p0, Vertex2d& p, int iq) const;

    public:
      Q42d();
      ~Q42d();

      std::string GetName() const {return "Q42d";}

      void BasicInit(const ParamFile* paramfile);
  };

/**********************************************************/

}

#endif
