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


#ifndef __hierarchicalmesh3d_h
#define __hierarchicalmesh3d_h

#include  "vertex.h" 
#include  "boundaryquad.h" 
#include  "hex.h" 
#include  "hexlawandorder.h" 
#include  "boundaryfunction.h"
#include  "hangcontainer3d.h" 
#include  "hierarchicalmesh.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class HierarchicalMesh3d : public HierarchicalMesh
{
  protected :
  
  /*  typedef  */

  typedef  std::vector<Vertex3d>       VertexVec3d;
  typedef  BoundaryCell<4>          BoundaryQuad;
  typedef  std::vector<Hex>            HexVec;
  typedef  std::vector<BoundaryQuad>   BQuadVec;
  typedef  HangList<2>              LineHangList;
  typedef  HangList<4>              QuadHangList;
  typedef  BoundaryFunction<3>      BoundaryFunction3d;
  typedef  std::map<int,fixarray<8,int> >  HexChilds;

  /*  Data  */

  CurvedShapes<3>    _curvedshapes;

  VertexVec3d        vertexs3d; 

  /* info fuer interpolation auf neues gitter */
  HexChilds          hexchildsofdeleted;
  HexVec             hexs;
  BQuadVec           Bquads;
  LineHangList       LineHang;
  QuadHangList       QuadHang;
  HexLawAndOrder     HexLaO;
  std::map<int,int>       hexofcurved;

  /*  Functionen  */
  int    Vater(const int i) const;
  IntVector Nachkommen(const int i) const;
  IntVector Geschwister(const int i) const;
  IntVector Kinder     (const int i) const;
  
  void post_refine3d();

  void delete_vertexs3d(const IntVector&);

  void new_edge_vertex3d(int, const EdgeVector&);
  void new_face_vertex3d(int, const FaceVector&);
  void new_vertex3d     (int, const fixarray<6,int>&);

  void check_mesh3d() const;

  std::pair<bool,tint>  check_inp(const std::string&);
  std::pair<int,int>  GetBoundaryInformation(int i) const;

  void init_quad           (BoundaryQuad&);

  void  build_neighbours() const;

  void prepare3d  (const IntVector&, const IntVector&, IntSet&, IntSet&);
  void new_hexs   (const HangContainer3d&, const IntVector&, 
		   const IntVector&, int, const IntSet&);
  void ghost_fill_neighbours2d();
  void ghost_fill_neighbours3d();
  void UpdateHangs(HangContainer3d& hangset,
		   const IntSet& cellref, 
		   const IntSet& cellcoarse);
  void FaceCoarse(HangContainer3d&, const IntSet&) const;
  void FaceRefine(HangContainer3d&, const IntSet&) const;
  void UpdateHangingEdges(HangContainer3d& hangset,
			  const IntSet& cellref, 
			  const IntSet& cellcoarse) const;
  void boundary_prepare3d(IntSet&, IntSet&, IntSet&, const HangContainer3d&);
  void new_boundary3d    (IntSet&, IntSet&,IntSet&);
  void new_vertexs3d     (HangContainer3d&, const IntVector&, const IntSet&);
  void basic_refine3d    (HangContainer3d&, const IntSet&, const IntSet&);
  void basic_fill_neighbours3d();
  void boundary_newton3d      (IntSet&);
  void inner_vertex_newton3d (const IntVector&, const IntSet&, const IntSet&);
  void update_boundary_data3d(const IntSet&);
  void new_bquads            (const IntVector&, const IntVector&, const IntSet&);
  void new_middle_vertex3d   (int,int);

  int  regular_grid3d_one  (IntSet&, IntVector&, const IntSet&, const IntSet&);
  int  regular_grid3d_one  (IntVector&, IntVector&, const IntSet&, const IntSet&);
  int  regular_grid3d_two  (IntVector&, const IntSet&);
  int  regular_grid3d_three_refine(IntSet&) const;
  int  regular_grid3d_three_coarse(IntSet&, IntSet&) const;

  void GetMinMaxLevels(IntVector& maxi, IntVector& mini,
		       const IntSet& CellRef) const;

  void init_edges3d();
  void LoadFathers3d(IntVector& v) const;

  void    _refine3d (IntSet&, IntSet&, const IntVector&, const IntVector&);
  void FillNeighbourFaces(const Hex& father, const FaceVector& Face,
			  int rneigh);
  void FillNeighbourFaces(int M, int S, const FaceVector& Face);
  void   InitHexOfCurved();
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

  HierarchicalMesh3d();
  HierarchicalMesh3d(const HierarchicalMesh3d& H);
  HierarchicalMesh3d& operator=(const HierarchicalMesh3d& H);
  HierarchicalMesh3d(const ParamFile* paramfile);
  ~HierarchicalMesh3d() {  GetCurvedShapes().clear();}

  std::string GetName() const {return "HierarchicalMesh3d";}

  /*  Zugriff  */

  int  dimension()            const { return 3;}

  int  nnodes   ()            const { return vertexs3d.size();}
  int  ncells   ()            const { return hexs.size();}
  int  nbquads  ()            const { return Bquads.size();}

  int  nodes_per_cell(int i)  const { return 8;}
  int  VtkType(int i) const { return 12;}

  const CurvedShapes<3>& GetCurvedShapes() const { return _curvedshapes;}
  CurvedShapes<3>& GetCurvedShapes() { return _curvedshapes;}

  const VertexVec3d& GetVertexVector() const {return vertexs3d; }
  VertexVec3d& GetVertexVector() {return vertexs3d; }

  const Vertex3d& vertex3d(int i)         const { return vertexs3d[i];}

  const Hex&   hex (int i)                const { return hexs[i];}
  const BoundaryQuad&  bquad(int i)       const { return Bquads[i];}

  int  vertex_of_cell (int i, int ii)      const { return hexs[i].vertex(ii);}
  int  vertex_of_bquad(int i, int ii)      const { return Bquads[i].vertex(ii);}
  int  face_of_hex    (int i, int ii)      const { return hexs[i].edge(ii); }
  int  level(int i)                        const { return hexs[i].level();}
  bool sleep(int i)                        const { return hexs[i].sleep();}

  int child(int i, int ii) const { return hexs[i].child(ii); }
  int nchilds(int i)       const { return hexs[i].nchilds(); }

  const HexLawAndOrder&     HexLawOrder()     const { return HexLaO; }
  const LineHangList&       linehang()        const { return LineHang;}
  const QuadHangList&       quadhanglist()    const { return QuadHang; }
  const BoundaryFunction3d* quad_shape(int i) const;

  const std::vector<BoundaryQuad>& quad_list() const { return Bquads; }

  const VertexVec3d&        vertex3d() const {return vertexs3d;}
  const HexVec&             hex     () const {return hexs;}
  const BQuadVec&           bquad   () const {return Bquads;}
  const QuadHangList&       quadhang() const {return QuadHang;}
  const std::map<int,int>&       GetHexOfCurved() const {return hexofcurved;}

  /*  Functionen  */

  void   write    (const std::string&) const;
  void   write_gup(const std::string&) const;
  void   write_gip(const std::string&) const;

  void   WriteAll(const std::string&) const;

  void   write_inp(const std::string&) const;
  void   read_inp (const std::string&);
  void   read_gup (const std::string&);
  void   read_gip (const std::string&);

  void   global_coarse3d();

  void   refine      (const IntVector&, const IntVector&);
  void   patch_refine(IntVector&, IntVector&);
  //  int    smooth_edges();
  void   FillAllBoundaryLines();

  pint  EdgeNeighbour(int i, int e) const;

  int  NodeOnFace(int e) const;
  fixarray<4,int> ChildrenOfFace(int e) const;

  void  GetVertexesOfFace(fixarray<4,int>&, int) const;
  void  GetVertexesOfFace(fixarray<5,int>&, int) const;
  void GetAwakePatchs(std::set<int>&) const;
  void GetAwakeCells(std::set<int>&) const;
  void ConstructQ2PatchMesh(IntVector& pm) const;
  IntVector ConstructQ4Patch(int c) const;
  std::set<int> GetColors() const;

  
  int nactivedescendants(int i)      const;
  IntVector GetVertices(int c) const;

  int GetBoundaryCellOfCurved(int iq) const
    {
      std::map<int,int>::const_iterator p = hexofcurved.find(iq);
      if( p!=hexofcurved.end() ) return p->second;
      return -1;
    }
  void Testing();
  int neighbour(int c, int le) const;
  int neighbour_neighbour(int c, int le) const;

  void AddShape(int col, BoundaryFunction<3>* f) {
    GetCurvedShapes().AddShape(col,f);
  }
};
}

/*---------------------------------------------------*/

#endif
