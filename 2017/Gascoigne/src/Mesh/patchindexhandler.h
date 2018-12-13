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


#ifndef  __PatchIndexHandler_h
#define  __PatchIndexHandler_h

#include  "gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class PatchIndexHandler
{
protected:

  bool                    haspatch,hasq4patch;
  nvector<IntVector>      indexofpatch,indexofq4patch;
  nvector<IntVector>      patch2cell,q4patch2cell;
  
  int dim;

public:
  PatchIndexHandler() : haspatch(false), hasq4patch(false) { }
  ~PatchIndexHandler() { }

  int&                       GetDim()          { return dim;}
  bool&                      GetHasPatch()     { return haspatch;}
  bool&                      GetHasQ4Patch()   { return hasq4patch;}
  nvector<IntVector>&        GetIndex()        { return indexofpatch;}
  const nvector<IntVector>&  GetIndex()const   { return indexofpatch;}
  nvector<IntVector>&        GetIndexQ4()      { return indexofq4patch;}
  const nvector<IntVector>&  GetIndexQ4()const { assert(hasq4patch); return indexofq4patch;}

  IntVector&          GetPatch2Cell(int i)
    { assert(i<patch2cell.size()); return patch2cell[i];}
  const IntVector&    GetPatch2Cell(int i) const
    { assert(i<patch2cell.size()); return patch2cell[i];}
  IntVector&          GetQ4Patch2Cell(int i)
    { assert(i<q4patch2cell.size()); return q4patch2cell[i];}
  const IntVector&    GetQ4Patch2Cell(int i) const
    { assert(hasq4patch && i<q4patch2cell.size()); return q4patch2cell[i];}
  
        nvector<IntVector>& GetAllPatch2Cell()         { return patch2cell; }
  const nvector<IntVector>& GetAllPatch2Cell() const   { return patch2cell; }
        nvector<IntVector>& GetAllQ4Patch2Cell()       { return q4patch2cell; }
  const nvector<IntVector>& GetAllQ4Patch2Cell() const { assert(hasq4patch); return q4patch2cell; }
  
  int npatches()    const { return indexofpatch.size();}
  int nq4patches()  const { return indexofq4patch.size();}
  bool HasPatch()   const { return haspatch;}
  bool HasQ4Patch() const { return hasq4patch;}
  int  Dim()        const { return dim; }
  
  const IntVector& IndicesOfPatch(int i)   const { return indexofpatch[i];}
  const IntVector& IndicesOfQ4Patch(int i) const { assert(hasq4patch); return indexofq4patch[i];}
  IntVector Q2IndicesOfQ4Patch(int i) const;
  IntVector CoarseIndices(int iq) const;
  IntVector CoarseIndicesQ4(int iq) const;

  int nodes_per_patch() const
    { 
      if (dim==2) return 9;
      return 27;
    }
  int nodes_per_q4patch() const
    { 
      if (dim==2) return 25;
      return 125;
    }
};
}

#endif
