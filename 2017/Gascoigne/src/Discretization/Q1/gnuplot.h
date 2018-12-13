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


#ifndef  __gnuplotdata_h
#define  __gnuplotdata_h

#include  <string>
#include "vertex.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class GnuplotData
{
protected:

  std::string   plane;
  Vertex3d pos;

public:

  GnuplotData() {};
  GnuplotData(const std::string& s, const Vertex3d& pos);
  GnuplotData(const GnuplotData& GP);
  
  void   SetName(std::string& filename) const;
  bool   TestVertex(const Vertex2d& v) const;
  bool   TestVertex(const Vertex3d& v) const;
  double SetVertex (const Vertex2d& v) const;
  double SetVertex (const Vertex3d& v) const;
};
}

#endif
