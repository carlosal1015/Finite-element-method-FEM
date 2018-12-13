/**
*
* Copyright (C) 2005, 2006 by the Gascoigne 3D authors
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


#ifndef __vtkvisu_h
#define __vtkvisu_h

#include  "visualization.h"

/*--------------------------------------------------*/

namespace Gascoigne
{

class VtkVisu : public Visualization
{
 public:

  VtkVisu(const MeshInterface& M, const std::string& name, int iter);
  ~VtkVisu();

  void WriteNodeData(const DoubleVector& eta);
  void WriteCellData(const GlobalVector& eta);
};

/*--------------------------------------------------*/
}

#endif
