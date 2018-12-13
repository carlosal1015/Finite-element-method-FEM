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


#include "regular_update.h"
#include <algorithm>


using namespace std;

/*----------------------------------------------*/

namespace Gascoigne
{
void regular_update(IntSet& hr, IntSet& hc, IntVector& vr, IntVector& vc)
{
  for (int i=0; i<vr.size(); i++)
    {
      hr.insert(vr[i]);
    }
  for (int i=0; i<vc.size(); i++)
    {
      hc.erase(vc[i]);
    }
  vr.insert(vr.end(),vc.begin(),vc.end());
  sort(vr.begin(),vr.end());
  int n = unique(vr.begin(),vr.end()) - vr.begin();
  vr.reserve(n); vr.resize(n);
}
}
