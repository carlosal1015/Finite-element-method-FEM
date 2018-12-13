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


#include  "columndiagstencil.h"

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne
{
void ColumnDiagStencil::memory(const SparseStructureInterface* SI)
{
  ColumnStencil::memory(SI);

  sdiag.reservesize(n());

  for(int i=0;i<n();i++)
    {
      for (int pos=sstart[i]; pos<sstart[i+1]; pos++)
        {
          if (scol[pos]==i)
            {
              sdiag[i] = pos;
              break;
            }
        }
    }
}

/*-------------------------------------------------------------*/

void ColumnDiagStencil::memory(int n, int nt)
{
  ColumnStencil::memory(n,nt);
  sdiag.reservesize(n);
}

}
