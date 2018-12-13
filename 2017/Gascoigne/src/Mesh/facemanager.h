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


#ifndef __facemanager_h
#define __facemanager_h

#include  "edge.h"
#include  "hexlawandorder.h" 
#include  "hangcontainer3d.h" 
#include  "gascoigne.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class FaceManager
{
  typedef fixarray<4,int>  FaceVector;

 protected:

  std::vector<Edge>&          edges;
  std::vector<Hex>&           hexs;
  const IntVector&            co2n;
        IntVector&            eo2n;

  IntVector                SwappedEdge;
  HexLawAndOrder           HexLaO;

  void Update      ();
  void InnerFaces  (const IntSet& CellRefList);
  void OuterFaces  (const HangContainer3d& hangset);
  void OldHangings (HangContainer3d& hangset3d, const IntSet& CellRefList);
  void SwappedFaces();
  void NeighbourTester() const;
  void FillNeighbourFaces(const Hex& M, const Hex& S, const FaceVector& Face);

 public:

  FaceManager(std::vector<Edge>&, std::vector<Hex>&, const IntVector& con, IntVector& eon);

  const Hex&  hex(int i)           const { return hexs[i];}
        Hex&  hex(int i)                 { return hexs[i];}

  fixarray<2,int> ChildrenOfFace(int e) const;

  bool EdgeIsHanging(int e) const;
  bool EdgeIsHanging(const Edge& e) const;

  void LoadFaceElimination(IntVector& edel, 
			   const IntSet& CellCoarseList,
			   const HangContainer3d& hangset) const;
  void Build( const IntSet& CellRefList, HangContainer3d&);
  void DeleteFaces();
  void InitFaces();
  void SortHangings();
  void Check(const HangContainer3d& hangset) const;
};
}

/*---------------------------------------------------*/

#endif
