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


#include  "quadlawandorder.h"

using namespace std;

/*----------------------------------------------------------------------*/

namespace Gascoigne
{
int QuadLawAndOrder::local_edge(const Quad& f, const EdgeVector& globaledge) const
{
  // globaledge muss sortiert sein !!!!!
  int li = f.global2local(globaledge[0]);
  // gives the local edge number of an edge

  LocVertexLocEdge::const_iterator lvle = lvlvle[li].begin();
  int   li2 = lvle->first;
  int   gi2 = f.vertex(li2); 

  if(find(globaledge.begin(),globaledge.end(),gi2)!=globaledge.end())
    {
      return lvle->second;
    }
  lvle++;
  return lvle->second;
}

/*----------------------------------------------------------------------*/

int QuadLawAndOrder::GlobalInnerEdge(int c, int i) const
{
  int ic = quads[c].child(i);
  int ie = InnerEdgeOfChild(i,0);
  return quads[ic].edge(ie);
}

/*----------------------------------------------------------------------*/

int QuadLawAndOrder::GlobalChildEdge(const EdgeVector& edge, int q, int j) const
{
  int ledge = local_edge_index(q,edge);
  int ic    = ChildsOfEdge(ledge,j);
  int child = quads[q].child(ic);
  int iedge = ChildEdge(ledge);
  
  return quads[child].edge(iedge);
}

/*---------------------------------------------------*/

pair<int,int> QuadLawAndOrder::GetChildEdges(EdgeVector& edge,
					     const EdgeVector& bigedge, 
					     int hanging, int bigquad, int i) const
{
  int bigeind = local_edge_index(bigquad,bigedge);
  int ic      = ChildsOfEdge(bigeind,i);
  int q       = quads[bigquad].child(ic);

  edge[0] = hanging; 
  edge[1] = bigedge[0];
  
  int iedge = local_edge_index(q,edge);
  
  if (iedge<0)
    {
      edge[1] = bigedge[1];
      iedge   = local_edge_index(q,edge);
    } 
  return make_pair(q,iedge);
}

/*----------------------------------------------------------------------*/

int  QuadLawAndOrder::global_index(const Quad& q, int i) const
{
  return quads[q.child(gc[i])].vertex(gv[i]);
}

/*----------------------------------------------------------------------*/

void QuadLawAndOrder::local_edge_index(EdgeVector& index, int edge) const
{
  index[0] = edge%4;
  index[1] = (edge+1)%4;
}

/*----------------------------------------------------------------------*/

void QuadLawAndOrder::globalvertices_of_edge(const Quad& q, EdgeVector& index, int edge) const
{
  local_edge_index(index, edge);

  index[0] = q.vertex(index[0]);
  index[1] = q.vertex(index[1]);
}

/*----------------------------------------------------------------------*/

void QuadLawAndOrder::fill_corner_vertex_in_childs(const Quad& f) const 
{
  quads[f.child(0)].vertex(0) = f.vertex(0);
  quads[f.child(1)].vertex(1) = f.vertex(1);
  quads[f.child(2)].vertex(2) = f.vertex(2);
  quads[f.child(3)].vertex(3) = f.vertex(3);
}

/*----------------------------------------------------------------------*/

void QuadLawAndOrder::fill_edge_vertex_in_childs(const Quad& f, int edge, 
						 int number) const
{
  int c0 = childs_edge[edge][0];
  int c1 = childs_edge[edge][1];
  quads[f.child(c0)].vertex(vice[edge][0]) = number;
  quads[f.child(c1)].vertex(vice[edge][1]) = number;
}

/*----------------------------------------------------------------------*/

void QuadLawAndOrder::fill_middle_vertex_in_childs(const Quad& f, 
						   int number) const  
{
  quads[f.child(0)].vertex(2) = number;
  quads[f.child(1)].vertex(3) = number;
  quads[f.child(2)].vertex(0) = number;
  quads[f.child(3)].vertex(1) = number;
}

/*----------------------------------------------------------------------*/

int QuadLawAndOrder::middle_vertex(const Quad& f) const
{
  return quads[f.child(0)].vertex(cell_midpoint(0));
}

/*----------------------------------------------------------------------*/

int QuadLawAndOrder::edge_vertex(const Quad& f, int edge) const
{
  if (!f.sleep()) return -1;
  int lic = childs_edge[edge][0];
  int liv = vice[edge][0];
  int gic = f.child(lic);

  return quads[gic].vertex(liv);
}

/*----------------------------------------------------------------------*/

void QuadLawAndOrder::childs_of_edge(QuadVector& child, const Quad& f,
				     int edge) const
{
  child[0] = f.child(childs_edge[edge][0]);
  child[1] = f.child(childs_edge[edge][1]);
}

/*----------------------------------------------------------------------*/

void QuadLawAndOrder::childs_of_global_edge(QuadVector& child, const Quad& f, 
					    const EdgeVector& globaledge) const
{
  int localedge = local_edge(f,globaledge);
  childs_of_edge(child,f,localedge);
}

/*----------------------------------------------------------------------*/

void QuadLawAndOrder::global_edge_unsorted(fixarray<2,int>& lineglob, const Quad& q,
				  int edge) const
{
  fixarray<2,int> lineloc;
  local_edge_index(lineloc,edge);
  q.vertex_loc2glob(lineglob,lineloc);
}

/*----------------------------------------------------------------------*/

void QuadLawAndOrder::globaledgechildren_of_father(vector<EdgeVector>& edges,
						   const Quad& f) const
{
  size_t n = 8;
  edges.resize(n);
  edges[0][0] = global_index(f,0); edges[0][1] = global_index(f,4);
  edges[1][0] = global_index(f,1); edges[1][1] = global_index(f,4);

  edges[2][0] = global_index(f,1); edges[2][1] = global_index(f,5);
  edges[3][0] = global_index(f,2); edges[3][1] = global_index(f,5);

  edges[4][0] = global_index(f,2); edges[4][1] = global_index(f,6);
  edges[5][0] = global_index(f,3); edges[5][1] = global_index(f,6);

  edges[6][0] = global_index(f,3); edges[6][1] = global_index(f,7);
  edges[7][0] = global_index(f,0); edges[7][1] = global_index(f,7);
      
  for(int i=0;i<n;i++) sort(edges[i].begin(),edges[i].end());
}

/*----------------------------------------------------------------------*/

int QuadLawAndOrder::local_edge_index(int q, const EdgeVector& edge) const
{
  for (int i=0; i<4; i++)
    {
      if ((edge[0]==quads[q].vertex(i)) && (edge[1]==quads[q].vertex((i+1)%4)))
	return i;
      if ((edge[1]==quads[q].vertex(i)) && (edge[0]==quads[q].vertex((i+1)%4)))
	return i;
    }
  return -1;
}

/*----------------------------------------------------------------------*/

QuadLawAndOrder::QuadLawAndOrder(vector<Quad>& q) : 
  quads(q),
  childs_edge(EdgeVector(2)), 
  vice(EdgeVector(2)), 
  child_point_cell(4), 
  child_point_vertex(4)
{
  vice[0][0] = 1; vice[0][1] = 0;
  vice[1][0] = 2; vice[1][1] = 1;
  vice[2][0] = 3; vice[2][1] = 2;
  vice[3][0] = 0; vice[3][1] = 3;
  
  childs_edge[0][0] = 0; childs_edge[0][1] = 1;
  childs_edge[1][0] = 1; childs_edge[1][1] = 2;
  
  childs_edge[2][0] = 2; childs_edge[2][1] = 3;
  childs_edge[3][0] = 3; childs_edge[3][1] = 0;
  
  child_point_cell[0] = 0; child_point_vertex[0] = 0;
  child_point_cell[1] = 1; child_point_vertex[1] = 1;
  child_point_cell[2] = 2; child_point_vertex[2] = 2;
  child_point_cell[3] = 3; child_point_vertex[3] = 3;

  gc[0] = 0;  gv[0] = 0;
  gc[1] = 1;  gv[1] = 1;
  gc[2] = 2;  gv[2] = 2;
  gc[3] = 3;  gv[3] = 3;
  gc[4] = 0;  gv[4] = 1;
  gc[5] = 1;  gv[5] = 2;
  gc[6] = 2;  gv[6] = 3;
  gc[7] = 3;  gv[7] = 0;
  gc[8] = 0;  gv[8] = 2;

  lvlvle.resize(4);
  lvlvle[0].insert(make_pair<int,int>(1,0));
  lvlvle[0].insert(make_pair<int,int>(3,3));
  
  lvlvle[1].insert(make_pair<int,int>(2,1));
  lvlvle[1].insert(make_pair<int,int>(0,0));
  
  lvlvle[2].insert(make_pair<int,int>(1,1));
  lvlvle[2].insert(make_pair<int,int>(3,2));
  
  lvlvle[3].insert(make_pair<int,int>(2,2));
  lvlvle[3].insert(make_pair<int,int>(0,3));

  ieoc[0][0] = 1; ieoc[0][1] = 2;
  ieoc[1][0] = 2; ieoc[1][1] = 3;
  ieoc[2][0] = 3; ieoc[2][1] = 0;
  ieoc[3][0] = 0; ieoc[3][1] = 1;

  oeoc[0][0] = 0; oeoc[0][1] = 3;
  oeoc[1][0] = 0; oeoc[1][1] = 1;
  oeoc[2][0] = 1; oeoc[2][1] = 2;
  oeoc[3][0] = 2; oeoc[3][1] = 3;

//   Pilot p;
//   p.check();
}
}
