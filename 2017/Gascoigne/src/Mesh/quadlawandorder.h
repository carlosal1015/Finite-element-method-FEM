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


#ifndef  __quadlawandorder_h
#define  __quadlawandorder_h

#include  "quad.h"
#include  <map>

  /*            2

    3__________6__________2
    |          |          |
    |          |          |
    |    3     |    2     |
    |          |          |
    |          |          |
3   7__________8__________5     1
    |          |          |
    |          |          |
    |    0     |    1     |
    |          |          |
    |          |          |
    0__________4__________1

               0
  */

/*---------------------------------------------------*/

namespace Gascoigne
{
class QuadLawAndOrder
{
 protected:

  typedef  std::map<int,int>               LocVertexLocEdge;
  typedef  std::vector<LocVertexLocEdge>   LocVertexLocVertexLocEdge;
  typedef fixarray<2,int>             EdgeVector;
  typedef fixarray<2,int>             QuadVector;

  // Daten fuer Suchen von kindern an hang

  LocVertexLocVertexLocEdge           lvlvle;

  // Referenzen fuer globale fkts

  std::vector<Quad>&          quads;

  fixarray<4,EdgeVector>   childs_edge, vice;
  fixarray<4,int>          child_point_cell, child_point_vertex;
  fixarray<9,int>          gc,gv;

  fixarray<4,fixarray<2,int> > ieoc, oeoc;

  int local_edge(const Quad& f, const EdgeVector& globaledge) const;

 public:

  QuadLawAndOrder(std::vector<Quad>& q);

  int  cell_midpoint    (int i)    const { return (i+2)%4; }
  int  global_index     (const Quad& q, int i) const;

  void fill_corner_vertex_in_childs (const Quad& f) const;
  void fill_edge_vertex_in_childs   (const Quad& f, int e, int i) const;
  void fill_middle_vertex_in_childs (const Quad& f, int i) const;

  //   edges

  int  ChildEdge        (int e)        const { return e; }
  int  ChildsOfEdge     (int e, int i) const { return childs_edge[e][i]; }
  int  InnerEdgeOfChild (int c, int i) const { return ieoc[c][i]; }
  int  OuterEdgeOfChild (int c, int i) const { return oeoc[c][i]; }
  int  GlobalInnerEdge(int c, int i) const;

  std::pair<int,int> GetChildEdges(EdgeVector& edge,const EdgeVector& bigedge, 
			      int hanging, int bigquad, int i) const;

  int  GlobalChildEdge (const EdgeVector& edge, int q, int j) const;
  void local_edge_index(EdgeVector& index, int edge) const;
  int  local_edge_index(int, const EdgeVector&) const;

  int middle_vertex (const Quad& f) const;
  int edge_vertex   (const Quad& f, int edge) const;

  /* for boundaries */
  void childs_of_edge(QuadVector& child, const Quad& f, int edge) const;

  /* for regular */
  void childs_of_global_edge(QuadVector& child, const Quad& f, 
			     const EdgeVector& globaledge) const;

  void globaledgechildren_of_father(std::vector<EdgeVector>& edges, 
				    const Quad& f) const;
  
  /* for hierarchicalmesh / find in liehanglist */
  void global_edge_unsorted(fixarray<2,int>& lineglob, const Quad& q, int edge) const;

  /* fuer mginterpolator */
  void globalvertices_of_edge(const Quad&, EdgeVector&, int) const;
};
}

#endif
