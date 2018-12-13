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


#include  "twinstencil.h"

using namespace std;

/*-------------------------------------------------------------*/
  
namespace Gascoigne
{
void TwinStencil::diagonalfirst()
{
  // has to be called if later on an ilu factorization is performed

  for(int i=0; i<n(); i++)
    {
      int first = sstart[i];
      int found = -1;
      for (int pos=first; pos<sstart[i+1]; pos++)
        {
          if (scol[pos]==i)
            {
              found=pos;
              break;
            }
        }
            if (found==-1)
        {
          cout << "UnstructuredStencil::diagonal not found " << i << endl;
          abort();
        }
            for (int pos=found; pos>first; pos--)
        {
          swap(scol[pos],scol[pos-1]);
        }
    }
}

/*-------------------------------------------------------------*/
  
int TwinStencil::half(int i) const
{
  for (int pos=sstart[i]; pos<sstart[i+1]; pos++)
    {
      if (scol[pos]>i) return pos;
    }
  return sstart[i+1];
}

/*-------------------------------------------------------------*/

void TwinStencil::memory(int n, int nt)
{
  ColumnStencil::memory(n,nt);
}

/*-------------------------------------------------------------*/

void TwinStencil::memory(const SparseStructureInterface* SA)
{
  ColumnStencil::memory(SA);
  diagonalfirst();
}
}
