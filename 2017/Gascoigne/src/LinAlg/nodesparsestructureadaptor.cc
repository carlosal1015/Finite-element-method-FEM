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


#include  "nodesparsestructureadaptor.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
void NodeSparseStructureAdaptor::FillStencil(ColumnDiagStencil& S) const
{
  S.start(0) = 0;
  for(int i=0;i<SSP->n();i++)
    {
      int srowsize = _ncomp*SSP->rowsize(i);
      for(int c=0;c<_ncomp;c++)
        {
          int ii = index(i,c);
          int first = S.start(ii);
          // 	  S.start(ii+1) = first + srowsize;
          S.stop(ii) = first + srowsize;
          
          int id = 0;
          for(SparseStructure::const_iterator p=SSP->rowbegin(i);
              p!=SSP->rowend(i);p++)
            {
              for(int d=0;d<_ncomp;d++)
                {
                  S.col(first+id) = index(*p,d);
                  id++;
                }
            }
        }
    }
}

/*-----------------------------------------*/

IntVector NodeSparseStructureAdaptor::GetIndicesDirichlet(int inode, const vector<int>& cv) const
{
  IntVector indices(cv.size());
  for(int ic=0;ic<cv.size();ic++)
    {
      indices[ic] = index(inode,cv[ic]);
    }
  return indices;
}
}
