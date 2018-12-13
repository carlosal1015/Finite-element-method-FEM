/**
*
* Copyright (C) 2004, 2005, 2007 by the Gascoigne 3D authors
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


#ifndef  __GascoigneMesh3d_h
#define  __GascoigneMesh3d_h

#include  "gascoignemesh.h"
#include  "vertex.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class GascoigneMesh3d : public GascoigneMesh
{
protected:

  // basic
  std::vector<Vertex3d>   nx;

public:

  GascoigneMesh3d();
  ~GascoigneMesh3d() {}

  std::string GetName() const {return "GascoigneMesh3d";}

        std::vector<Vertex3d>& GetVertexVector()       {return nx;}
  const std::vector<Vertex3d>& GetVertexVector() const {return nx;}

  int  dimension() const {return 3;}
  int  nnodes()    const {return nx.size();}
  int  ncells()    const {return nc.size()/8;}
  int  nhanging()  const { return HangingHandler.GetStructure()->size()
			     + HangingHandler.GetStructureFace()->size(); }

  int  nodes_per_cell(int i)  const { return 8;}
  int  VtkType(int i) const { return 12;}

  const Vertex3d& vertex3d(int i) const { return nx[i];} 
  int  vertex_of_cell(int i, int ii) const { return nc[8*i+ii]; }

  IntVector  IndicesOfCell(int iq) const;
};
}

#endif
