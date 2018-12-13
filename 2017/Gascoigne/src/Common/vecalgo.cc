/**
*
* Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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


#include  "vecalgo.h"
#include  "gascoignemath.h"
#include  <algorithm>
#include  "giota.h"

using namespace std;

namespace Gascoigne
{

/*************************************************************/

void transfer(int n, vector<int>& tr, const set<int>& del)
{
  tr.resize(n, -1);

  int count = 0;
  for(int i = 0; i < n; ++i)
  {
    if(del.find(i) == del.end())
    {
      tr[i] = count++;
    }
  }
}

/*************************************************************/

void transfer(int n, vector<int>& tr, vector<int>& del)
{
  tr.resize(n, -1);

  if (del.size() == 0)
  {
    iota(tr.begin(), tr.end(), 0);
    return;
  }

  sort(del.begin(),del.end());

  int count = 0;
  int pos   = 0;

  for(int i = 0; i < n; ++i)
  {
    while ((pos < del.size()) && (del[pos] < i))
    {
      pos++;
    }

    if ((pos == del.size()) || (del[pos] > i))
    {
      tr[i] = count++;
    }
  }
}
}
