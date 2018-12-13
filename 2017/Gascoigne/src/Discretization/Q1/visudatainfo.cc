/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#include  "visudatainfo.h"
#include  "compose_name.h"

using namespace std;

/*-------------------------------------------------------------------------*/

namespace Gascoigne
{
VisuDataInfo& VisuDataInfo::operator=(const VisuDataInfo& V)
{
  scalars = V.Scalars();
  vectors = V.Vectors();
  return *this;
}

/*-------------------------------------------------------------------------*/

bool VisuDataInfo::operator!=(const VisuDataInfo& V) const
{
  return (scalars!=V.Scalars())||(vectors!=V.Vectors());
}

/*-------------------------------------------------------------------------*/

VisuDataInfo::VisuDataInfo(const VisuData& D, string def)
{
  for(int c=0;c<D.visucomp();c++)
    {
      string name(def);
      compose_name_without_dot(name,c);
      AddScalar(c,name,c);
    }
}

/*-------------------------------------------------------------------------*/

void VisuDataInfo::AddScalars(int ncomp, string def)
{
  for(int c=0;c<ncomp;c++)
    {
      string name(def);
      compose_name_without_dot(name,c);
      AddScalar(c,name,c);
    }
}
}
