/**
*
* Copyright (C) 2004, 2007, 2009 by the Gascoigne 3D authors
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


#ifndef __boundaryindexhandler_h
#define __boundaryindexhandler_h

#include <map>
#include <set>
#include  "gascoigne.h"

/*--------------------------------------------------------------*/

namespace Gascoigne
{
class BoundaryIndexHandler
{
 protected:

  typedef std::map<int,IntVector> VecMap;

  IntSet  AllColors;
  VecMap  verteces, cells, localci, patches, localpi;
  
  std::map<int,std::map<int,int> > _PeriodicPairs;

 public:

  void CopySetToVector(const std::vector<IntSet>&,
		       const IntVector&, VecMap&) const;

  void clear();

  const IntSet& GetColors() const {return AllColors;}
  const VecMap& GetVertex() const {return verteces;}
  const VecMap& GetCell()   const {return cells;}
  const VecMap& GetLocal() const  {return localci;}
  const VecMap& GetPatch()   const {return patches;}
  const VecMap& GetLocalPatch() const  {return localpi;}

  IntSet& GetColors()  {return AllColors;}
  VecMap& GetVertex()  {return verteces;}
  VecMap& GetCell()    {return cells;}
  VecMap& GetLocal()   {return localci;}
  VecMap& GetPatch()   {return patches;}
  VecMap& GetLocalPatch() {return localpi;}

  const IntVector& Verteces(int col) const;
  const IntVector& Cells   (int col) const;
  const IntVector& Localind(int col) const;
  const IntVector& Patches   (int col) const;
  const IntVector& LocalPatchind(int col) const;

  void SetPeriodicPairs(std::map<int,std::map<int,int> > mm_PeriodicPairs);
  const std::map<int,std::map<int,int> > GetPeriodicPairs() const;

  void Equal(const IntSet& col, const VecMap& v, const VecMap& c, const VecMap& l);
  void check() const;
  friend std::ostream& operator<<(std::ostream &s, const BoundaryIndexHandler& A);
};
}

/*--------------------------------------------------------------*/

#endif
