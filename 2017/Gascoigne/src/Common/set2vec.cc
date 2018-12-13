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


#include  "set2vec.h"

using namespace std;

namespace Gascoigne
{

/*---------------------------------------------------*/

void Set2Vec(vector<int>& v, const set<int>& h)
{
  v.resize(h.size());
  int j = 0;
  for (set<int>::const_iterator p=h.begin();
       p!=h.end(); p++)
    {
      v[j++] = *p;
    }
}

/*---------------------------------------------------*/

void Vec2Set(set<int>& h, const vector<int>& v)
{
  h.clear();
  for (int i=0; i<v.size(); i++)
    {
      h.insert(v[i]);
    }
}
}
