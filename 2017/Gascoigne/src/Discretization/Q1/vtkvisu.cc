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


#include "vtkvisu.h"
#include "visudatanvector.h"
#include "visudatacompvector.h"

namespace Gascoigne
{
/*-------------------------------------------------*/

VtkVisu::~VtkVisu() {}

/*-------------------------------------------------*/

VtkVisu::VtkVisu(const MeshInterface& M, const std::string& name, int iter) : Visualization()
{
  format("vtk");
  set_name(name);
  step(iter);

  SetMesh(&M);
}

/*-------------------------------------------------*/

void VtkVisu::WriteNodeData(const DoubleVector& eta)
{
  VisuDataInfo     VDI(1);
  VisuDataNVector  VD(eta);
  
  SetPointData(&VD);
  SetPointDataInfo(&VDI);

  write();
}

/*-------------------------------------------------*/

void VtkVisu::WriteCellData(const GlobalVector& eta)
{
  int ncomp = eta.ncomp();

  VisuDataInfo     VDI(1);
  VDI.Clear();
  //  VisuDataNVector  VD(eta);
  VisuDataCompVector VD(eta);

  VDI.AddScalars(ncomp);

  SetCellData(&VD);
  SetCellDataInfo(&VDI);

  write();
}

/*-------------------------------------------------*/
}
