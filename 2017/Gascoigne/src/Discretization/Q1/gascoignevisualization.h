/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#ifndef  __GascoigneVisualization_h
#define  __GascoigneVisualization_h

#include  "gascoigne.h"
#include  "visualization.h"
#include  "visudatacompvector.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class ComponentInformation;

class GascoigneVisualization : public Visualization
{
protected:

  const GlobalVector* _v;

  VisuDataInfo        VDI;
  VisuDataCompVector  VD;

  void AddVector(const GlobalVector* v);
  void AddVector(const ComponentInformation* CI, const GlobalVector* v);

public:

  GascoigneVisualization() : Visualization(), _v(NULL) {}
  ~GascoigneVisualization() {}

  void AddPointVector(const ComponentInformation* CI, const GlobalVector* v);
  void AddPointVector(const GlobalVector* v);
  void AddCellVector(const ComponentInformation* CI, const GlobalVector* v);
  void AddCellVector(const GlobalVector* v);
};
}

#endif
