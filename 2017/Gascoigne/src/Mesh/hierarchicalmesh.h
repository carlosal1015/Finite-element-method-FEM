/**
*
* Copyright (C) 2004, 2005, 2006, 2011 by the Gascoigne 3D authors
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


#ifndef __hierarchicalmesh_h
#define __hierarchicalmesh_h

#include  <utility>

#include  "triple.h"
#include  <string>
#include  "edge.h" 
#include  "hanglist.h" 
#include  "curvedshapes.h"
#include  "meshinterface.h"
#include  "vertex.h"
#include  "compvector.h"
#include  "gascoigne.h"
#include  "paramfile.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class HierarchicalMesh : public MeshInterface
{
 protected :

  /*  typedef  */

  typedef  std::pair<int,int>      pint;
  typedef  triple<int,int,int>     tint;
  typedef  fixarray<2,int>         EdgeVector;
  typedef  fixarray<4,int>         FaceVector;
  typedef  std::vector<Edge>       EdgeVec;
  typedef  IntSet::iterator        IntSetIt;   
  typedef  IntSet::const_iterator  IntSetCIt;
 
  /*  Data  */

  int       mnlevels, pdepth, etapatcher;
  IntVector    vo2n, eo2n, co2n;
  EdgeVec   edges;

  mutable int _i_showoutput;
  
  void  update_edges(IntVector&);
  virtual int   FindPatchDepth() const=0;
  virtual void  FillVertexLevels(IntVector& dst) const=0;
  virtual void  RefineCoarseNodes(IntSet& dst, const IntVector& refnodes,
				  const IntVector& vertexlevel) const=0;
  virtual void  VertexToCells(IntVector& dst, const IntSet& src, 
			      const IntVector& vertexlevel) const=0;
  virtual void VertexToCellsCoarsening(IntVector& dst, const IntSet& src, 
			       const IntVector& vertexlevel) const=0;

  public:

  virtual ~HierarchicalMesh();
  
  int  withfaces;

  HierarchicalMesh();
  HierarchicalMesh(const HierarchicalMesh&);
  HierarchicalMesh& operator=(const HierarchicalMesh&);

  /*  Zugriff  */

  virtual int  nnodes ()        const =0;
  int  nlevels()                const { return 1+mnlevels;}
  int  nedges ()                const { return edges.size();}

  const IntVector*  Vertexo2n()    const { return &vo2n;}
  const IntVector*  Edgeo2n  ()    const { return &eo2n;}
  const IntVector*  Cello2n  ()    const { return &co2n;}

  int Vertexo2n(int i)          const { assert(i<vo2n.size());return vo2n[i];}
  int Edgeo2n  (int i)          const { assert(i<eo2n.size());return eo2n[i];}
  int Cello2n  (int i)          const { assert(i<co2n.size());return co2n[i];}

  const Edge&    edge(int i)    const {assert(i<edges.size());return edges[i];}
  const EdgeVec& edge()         const { return edges;}

  virtual int child(int i, int ii) const =0;
  virtual int nchilds(int i)       const =0;

  virtual int  level(int i)   const =0;
  virtual bool sleep(int i)   const =0;

  virtual int    Vater(const int i) const
    {
      std::cerr << "\"HierarchicalMesh::Vater\" not written!" << std::endl;
      abort();
    }
  virtual IntVector Nachkommen(const int i) const
    {
      std::cerr << "\"HierarchicalMesh::Nachkommen\" not written!" << std::endl;
      abort();
    }
  virtual IntVector Geschwister(const int i) const
    {
      std::cerr << "\"HierarchicalMesh::Geschwister\" not written!" << std::endl;
      abort();
    }
  virtual IntVector Kinder(const int i) const
    {
      std::cerr << "\"HierarchicalMesh::Kinder\" not written!" << std::endl;
      abort();
    }
    
  void SetParameters(std::string gridname, int patchdepth, int epatcher);
  void ReadFile(const std::string& gridname);
  void BasicInit(const ParamFile* pf);
  void global_refine  (int k);
  void global_patch_coarsen (int k);
  void random_refine  (double, int k=1);
  void random_patch_refine  (double, int k=1);
  void random_patch_coarsen (double, int k=0);
  void random_double_patch_refine  (double, int k=1);
  void clear_transfer_lists();
  virtual void   write_gip(const std::string&) const =0;
  virtual void   write_gup(const std::string&) const =0;
  virtual void   write_inp(const std::string&) const =0;

  virtual int dimension()       const { return 0;}
  virtual int ncells ()         const =0;
  int patchdepth()              const {return pdepth;}
  virtual int nactivedescendants(int i)      const=0;
  virtual IntVector GetVertices(int c) const=0;
  
  bool CellIsCurved(int iq) const  { return GetBoundaryCellOfCurved(iq)!=-1;}

  virtual std::set<int> GetColors() const=0;

  virtual void read_inp (const std::string&)=0;
  virtual void read_gup (const std::string&)=0;
  virtual void read_gip (const std::string&)=0;
  virtual void refine(const IntVector&, const IntVector&)=0;
  virtual void patch_refine(IntVector&, IntVector&)=0;
  virtual void vertex_patch_refine(IntVector& ref, IntVector& coarse);
  virtual void vertex_patch_refine(IntVector&);
  virtual void GetAwakePatchs(std::set<int>&) const =0;
  virtual void GetAwakeCells(std::set<int>&) const =0;
  virtual void ConstructQ2PatchMesh(IntVector& pm) const=0;
  virtual IntVector ConstructQ4Patch(int c) const =0;
  virtual std::set<int> CellNeighbours(int i) const 
    { std::cerr << "no CellNeighbours"; abort(); return std::set<int>();}

  virtual int  GetBoundaryCellOfCurved(int iq) const { return -1;}

  virtual void AddShape(int col, BoundaryFunction<2>* f)
  {
    std::cerr << "\"HierarchicalMesh::AddShape\" not written!" << std::endl;
    abort();
  }
  virtual void AddShape(int col, BoundaryFunction<3>* f)
  {
    std::cerr << "\"HierarchicalMesh::AddShape\" not written!" << std::endl;
    abort();
  }

  void ShowOutput(int i) const { _i_showoutput=i; } 

};
}

/*---------------------------------------------------*/

#endif
