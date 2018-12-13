/**
*
* Copyright (C) 2004, 2006, 2007 by the Gascoigne 3D authors
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


#ifndef  __levelmesh3d_h
#define  __levelmesh3d_h

#include  "hierarchicalmesh3d.h"
#include  "index.h"
#include  "boundaryindexhandler.h"
#include  "patchindexhandler.h"
#include  "gascoigne.h"

/*--------------------------------------------------------------*/

namespace Gascoigne
{
class LevelMesh3d : public Index
{
 protected:

  typedef std::map<int,fixarray<3,int> >  QuadraticHNStructure3;
  typedef std::map<int,fixarray<9,int> >  QuadraticHNStructure9;

  typedef std::map<int,fixarray<6,int> >   QuarticHNStructure5;
  typedef std::map<int,fixarray<26,int> >  QuarticHNStructure25;

  const HierarchicalMesh3d* HMP;

  void check_leveljump() const;
  void fill_opis(IntSet& dst, IntSet& oldquads) const;
  void fill_enkel(IntSet& dst, const Hex& Q) const;
  void fill_childs(IntSet& dst, const Hex& Q) const;
  bool EnkelUniform(const Hex& Q) const;
  bool BuildFathers(std::set<int>&  Vaeter) const;
  int  hexindex  (int i) const { return Index::Hexl2g(i);}
  void InitCells(int n);
  void InitNodes(int n);
  void InitEdges(int n);

  int refine_level(int n) const;
  void ConstructNodesOnFaceQ4(fixarray<81,int>& nodesonface,int vater,int ni) const;
  void InsertHangingFacesQ4(QuarticHNStructure25& hnq4face,const fixarray<81,int>& nodesonface) const;
  void InsertHangingEdgesQ4(QuarticHNStructure5&  hnq4, const fixarray<81,int>& nodesonface) const;
  void ConstructNodesOnFace(fixarray<25,int>& nodesonface,int vater,int ni) const;
  void InsertHangingFaceQ4 (QuarticHNStructure25& hnq4face,const fixarray<81,int>& nodesonface,
			    int n1,int n2,int n3,int n4,const fixarray<25,int>& I) const;
  void InsertHangingEdgeQ4 (QuarticHNStructure5&   hnq4,const fixarray<81,int>& nof,
			    int n1,int n2,int n3,int n4,int i1,int i2,int i3,int i4,int i5) const;

 public:
   
  LevelMesh3d(const HierarchicalMesh* hmp);
  ~LevelMesh3d();

  int ncells   ()         const  { return Index::HexSize(); }
  int dimension()         const  { return HMP->dimension();}

  int vertex_of_cell(int i, int ii) const 
    { return Index::Vertexg2l(HMP->vertex_of_cell(hexindex(i),ii)); }

  const Vertex3d& vertex3d(int i) const 
    { return HMP->vertex3d(Index::Vertexl2g(i)); }

  /*----- Functions -----*/

  const Hex&   hex  (int i) const { return HMP->hex(hexindex(i));}

  void BasicInit(const IntSet& n, const IntSet& o);

  /*----- Functions fuer Patch-----*/

  void construct_lists(IntSet& newquads, IntSet& oldquads) const;

  void ConstructHangingStructureQuadratic(QuadraticHNStructure3& hnq2,
					  QuadraticHNStructure9& hnq2face) const;
  void ConstructHangingStructureQuartic(QuarticHNStructure5& hnq4,
					QuarticHNStructure25& hnq4face) const;

  void InitBoundaryHandler(BoundaryIndexHandler& BI,const PatchIndexHandler& PIH) const;
  bool ConstructCellIndOfPatch(IntVector& dstc) const;
  void ConstructIndOfPatch(nvector<IntVector>& dstv) const;
};
}

/*---------------------------------------------------*/

#endif
