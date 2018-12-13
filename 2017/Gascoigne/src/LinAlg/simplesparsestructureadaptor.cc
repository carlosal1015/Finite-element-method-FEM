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


#include  "simplesparsestructureadaptor.h"

/*-----------------------------------------*/

namespace Gascoigne
{
void SimpleSparseStructureAdaptor::FillStencil(ColumnDiagStencil& S) const
{
  S.start(0) = 0;
  for(int i=0;i<SSP->n();i++)
    {
      int srowsize = SSP->rowsize(i);
      int first = S.start(i);
      S.stop(i) = first + srowsize;
	  
      int id = 0;
      for(SparseStructure::const_iterator p=SSP->rowbegin(i);
	  p!=SSP->rowend(i);p++)
	{
	  S.col(first+id++) = *p;
	}
    }
}
}
