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


#include  "columnstencil.h"
#include  "sparsestructure.h"


using namespace std;

/*-------------------------------------------------------------*/
  
namespace Gascoigne
{
ostream& operator<<(ostream &s, const ColumnStencil& A)
{
  s << "start:\n"<< A.start() << endl;
  s << "col:\n"<< A.col() << endl;
  return s;
}

/*-------------------------------------------------------------*/

void ColumnStencil::memory(int n, int nt)
{
  scol  .reservesize(nt);
  sstart.reservesize(n+1);
}

/*-------------------------------------------------------------*/

void ColumnStencil::memory(const SparseStructureInterface* SI)
{
  const SparseStructure* SS = dynamic_cast<const SparseStructure*>(SI);
  assert(SS);

  memory(SS->n(),SS->ntotal());

  sstart[0] = 0;
  for(int i=0;i<SS->n();i++)
    {
      int first = sstart[i];
      sstart[i+1] = first + SS->rowsize(i);
      int id = 0;
      for(set<int>::const_iterator p=SS->rowbegin(i);
          p!=SS->rowend(i);p++)
        {
          scol[first+id] = *p;
          id++;
        }
    }
}

/*-------------------------------------------------------------*/

std::ostream& ColumnStencil::Write(std::ostream& os) const
{
  os << n() << "\t" << nentries() << "\n\n"<< sstart << "\n\n";
  for(int i=0;i<n();i++)
    {
      for(int pos=start(i);pos<stop(i);pos++)
	{
	  os << col(pos) << " ";
	}
      os << std::endl;
    }
  return os;
}

/*-------------------------------------------------------------*/

}
