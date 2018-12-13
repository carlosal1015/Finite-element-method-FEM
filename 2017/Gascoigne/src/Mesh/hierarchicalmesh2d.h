/**
*
* Copyright (C) 2004, 2005, 2006, 2007 by the Gascoigne 3D authors
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


#ifndef __hierarchicalmesh2d_h
#define __hierarchicalmesh2d_h

#include  "vertex.h" 
#include  "boundaryline.h" 
#include  "quad.h" 
#include  "quadlawandorder.h" 
#include  "boundaryfunction.h"
#include  "hangcontainer2d.h"
#include  "hierarchicalmesh.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class HierarchicalMesh2d : public HierarchicalMesh
{
  protected :
  
  typedef  std::vector<Vertex2d>       VertexVec2d;
  typedef  std::vector<Quad>           QuadVec;
  typedef  BoundaryFunction<2>      BoundaryFunction2d;
  typedef  BoundaryCell<2>          BoundaryLine;
  typedef  std::vector<BoundaryLine>   BLineVec;
  typedef  HangList<2>              LineHangList;

  typedef triple<int,int,int>          tint;
  typedef std::map<fixarray<2,int>,HierarchicalMesh2d::BoundaryLine> HangBLList;

  /*  Data  */

  CurvedShapes<2>    _curvedshapes;

  VertexVec2d        vertexs2d; 

  QuadVec            quads;
  BLineVec           Blines;
  LineHangList       LineHang;
  QuadLawAndOrder    QuadLaO;
  std::map<int,int>       quadofcurved;

  /*  Functionen  */

  void post_refine2d();

  void delete_vertexs2d(const IntVector&);

  void new_edge_vertex2d(int, const EdgeVector&);
  void new_face_vertex2d(int, const FaceVector&);

  void  check_mesh2d() const;
  void  prepare2d  (const IntVector&, const IntVector&, IntSet&, IntSet&);
  std::pair<bool,tint>  check_inp(const std::string&);
  void  ghost2d    (HangContainer2d&, const IntSet&, const IntSet&);
  void  ghostglobalcoarse(HangContainer2d&, const IntSet&);
  void  ghost_fill_neighbours2d();
  void  basic_fill_neighbours2d();
  void  new_vertexs2d(HangContainer2d&, const IntVector&, const IntSet&);
  void  new_quads(const HangContainer2d&, const IntVector&, const IntVector&, 
		  int, const IntSet&);

  void  change_hangs2d  (const IntVector&, const IntVector&);
  void  change_vertexs2d(const IntVector&);
  void  change_quads2d  (const IntVector&, const IntVector&);
  void  boundary_prepare2d(IntSet&, IntSet&, IntSet&, const HangContainer2d&);
  void  new_boundary2d    (IntSet&, IntSet&, IntSet&);

  void  basic_refine2d(HangContainer2d&, const IntSet&, const IntSet&);

  void init_line           (BoundaryLine&);
  void new_lines           (const IntVector&, const IntVector&, const IntSet&);
  void boundary_newton2d     ();
  void inner_vertex_newton2d (const IntVector&, const IntSet&);
  void update_boundary_data2d(const IntSet&);

  int   regular_grid2d_one  (IntSet&, IntVector& , IntSet&, IntSet&) const;
  int   regular_grid2d_two  (IntSet&, IntSet&) const;
  int   regular_grid2d_three(IntSet&, IntSet&) const;
  int   regular_grid2d_three_refine(IntSet&) const;
  int   regular_grid2d_three_coarse(IntSet&, IntSet&) const;

  void GetMinMaxLevels(IntVector& maxi, IntVector& mini,
		       const IntSet& CellRef) const;
  void  init_edges2d();

  void LoadFathers  (IntVector& v) const;

  void    _refine2d (IntSet&, IntSet&, const IntVector&, const IntVector&);
  void   InitQuadOfCurved();
  int   FindPatchDepth() const;
  void  FillVertexLevels(IntVector& dst) const;
  void  RefineCoarseNodes(IntSet& dst, const IntVector& refnodes,
			  const IntVector& vertexlevel) const;
  void  VertexToCells(IntVector& dst, const IntSet& src, 
		      const IntVector& vertexlevel) const;
  void VertexToCellsCoarsening(IntVector& dst, const IntSet& src, 
			       const IntVector& vertexlevel) const;
  void recursive_childs(int q, IntVector& ref, int d) const;


  public:
  
  HierarchicalMesh2d();
  HierarchicalMesh2d(const HierarchicalMesh2d& H);
  HierarchicalMesh2d& operator=(const HierarchicalMesh2d& H);
  HierarchicalMesh2d(const ParamFile* paramfile);
  ~HierarchicalMesh2d()   {  GetCurvedShapes().clear();}

  std::string GetName() const {return "HierarchicalMesh2d";}

  /*  Zugriff  */

  int  dimension()           const { return 2;}
  int  nnodes   ()           const { return vertexs2d.size();}
  int  ncells   ()           const { return quads.size();}
  int  nblines  ()           const { return Blines.size(); }
  int  nodes_per_cell(int i)   const { return 4;}
  int VtkType(int i) const { return 9;}

  const CurvedShapes<2>& GetCurvedShapes() const { return _curvedshapes;}
  CurvedShapes<2>& GetCurvedShapes() { return _curvedshapes;}

  const VertexVec2d& GetVertexVector() const {return vertexs2d; }
  VertexVec2d& GetVertexVector() {return vertexs2d; }
  
  const Vertex2d& vertex2d(int i)         const { return vertexs2d[i];}

  const Quad&  quad(int i)                const { return quads[i];}
  const BoundaryLine&  bline(int i)       const { return Blines[i];}
  std::pair<int,int>  GetBoundaryInformation(int i) const;

  int  vertex_of_cell(int i, int ii)      const { return quads[i].vertex(ii);}
  int  vertex_of_bline(int i, int ii)     const { return Blines[i].vertex(ii);}
  int  edge_of_quad  (int i, int ii)      const { return quads[i].edge(ii); }
  int  level(int i)                       const { return quads[i].level();}
  bool sleep(int i)                       const { return quads[i].sleep();}

  int child(int i, int ii) const { return quads[i].child(ii); }
  int nchilds(int i)       const { return quads[i].nchilds(); }

  int   QuadNeighbour(const Quad&, int)    const;

  const QuadLawAndOrder&    QuadLawOrder()    const { return QuadLaO; }
  const LineHangList&       linehanglist()    const { return LineHang; }
  const BoundaryFunction2d* line_shape(int i) const;

  const std::vector<BoundaryLine>& line_list() const { return Blines; }

  const VertexVec2d&        vertex2d() const {return vertexs2d;}
  const QuadVec&            quad    () const {return quads;}
  const BLineVec&           bline   () const {return Blines;}
  const LineHangList&       linehang() const {return LineHang;}
  const std::map<int,int>&       GetQuadOfCurved() const {return quadofcurved;}

  /*  Functionen  */

  int    Vater(const int i) const;
  IntVector Nachkommen(const int i) const;
  IntVector Geschwister(const int i) const;
  IntVector Kinder     (const int i) const;
  int nactivedescendants(int i)      const;
  IntVector GetVertices(int c) const;
  
  void   write    (const std::string&) const;
  void   write_gup(const std::string&) const;
  void   write_gip(const std::string&) const;

  void   WriteAll(const std::string&) const;

  void   write_inp(const std::string&) const;
  void   read_inp (const std::string&);
  void   write_vtk(const std::string&) const;

  void   read_gup (const std::string&);
  void   read_gip (const std::string&);

  void   global_coarse();

  void   refine      (const IntVector&, const IntVector&);
  void   patch_refine  (IntVector&, IntVector&);
  int    smooth_edges();
  void   FillAllBoundaryLines();

  pint  EdgeNeighbour(int i, int e) const;
  void  VertexNeighbours2d(std::set<int>&, int i) const;

  int  NodeOnEdge(int e) const;
  fixarray<2,int> ChildrenOfEdge(int e) const;

  void  GetVertexesOfEdge(fixarray<3,int>&, int) const;
  void  GetVertexesOfEdge(fixarray<2,int>&, int) const;
  void GetAwakePatchs(std::set<int>&) const;
  void GetAwakeCells(std::set<int>&) const;
  void ConstructQ2PatchMesh(IntVector& pm) const;
  IntVector ConstructQ4Patch(int c) const;
  std::set<int> GetColors() const;
  int GetBoundaryCellOfCurved(int iq) const
    {
      std::map<int,int>::const_iterator p = quadofcurved.find(iq);
      if( p!=quadofcurved.end() ) return p->second;
      return -1;
    }

  std::set<int>   CellNeighbours(int i)    const;

  int neighbour(int c, int le) const;
  void FillVolumes(DoubleVector& vol) const;
  void writeq2(const IntVector &a,const std::vector<int> & b,int np) const;

  void AddShape(int col, BoundaryFunction<2>* f) {
    GetCurvedShapes().AddShape(col,f);
  }
};
}

/*---------------------------------------------------*/

#endif
