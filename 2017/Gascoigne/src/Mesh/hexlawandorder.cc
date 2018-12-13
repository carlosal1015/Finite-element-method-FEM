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


#include  "hexlawandorder.h"
#include  <algorithm>
#include  "edgearray.h"

using namespace std;

/*----------------------------------------------------------------------*/

namespace Gascoigne
{
HexLawAndOrder::HexLawAndOrder(vector<Hex>& h) : 
  hexs(h)
{
  // gibt die lokalen vertices einer face zurueck
  lvf[0][0] = 0; lvf[0][1] = 1; lvf[0][2] = 2;  lvf[0][3] = 3;
  lvf[1][0] = 1; lvf[1][1] = 5; lvf[1][2] = 6;  lvf[1][3] = 2;
  lvf[2][0] = 2; lvf[2][1] = 6; lvf[2][2] = 7;  lvf[2][3] = 3;
  lvf[3][0] = 3; lvf[3][1] = 7; lvf[3][2] = 4;  lvf[3][3] = 0;
  lvf[4][0] = 0; lvf[4][1] = 4; lvf[4][2] = 5;  lvf[4][3] = 1;
  lvf[5][0] = 4; lvf[5][1] = 7; lvf[5][2] = 6;  lvf[5][3] = 5;

  // gibt die lokalen vertices einer edge zurueck
  lve[0][0]  = 0; lve[0][1]  = 1;
  lve[1][0]  = 1; lve[1][1]  = 2;
  lve[2][0]  = 2; lve[2][1]  = 3;
  lve[3][0]  = 3; lve[3][1]  = 0;
  lve[4][0]  = 4; lve[4][1]  = 5;
  lve[5][0]  = 5; lve[5][1]  = 6;
  lve[6][0]  = 6; lve[6][1]  = 7;
  lve[7][0]  = 7; lve[7][1]  = 4;
  lve[8][0]  = 0; lve[8][1]  = 4;
  lve[9][0]  = 1; lve[9][1]  = 5;
  lve[10][0] = 2; lve[10][1] = 6;
  lve[11][0] = 3; lve[11][1] = 7;

  // gibt die lokalen childs einer face zurueck
  childs_face[0][0] = 0; childs_face[0][1] = 1;
  childs_face[0][2] = 2; childs_face[0][3] = 3;
  childs_face[1][0] = 1; childs_face[1][1] = 5;
  childs_face[1][2] = 6; childs_face[1][3] = 2;
  childs_face[2][0] = 3; childs_face[2][1] = 2;
  childs_face[2][2] = 6; childs_face[2][3] = 7;
  childs_face[3][0] = 4; childs_face[3][1] = 0;
  childs_face[3][2] = 3; childs_face[3][3] = 7;
  childs_face[4][0] = 4; childs_face[4][1] = 5;
  childs_face[4][2] = 1; childs_face[4][3] = 0;
  childs_face[5][0] = 5; childs_face[5][1] = 4;
  childs_face[5][2] = 7; childs_face[5][3] = 6;

  // gibt die lokalen childs an einer edge zurueck
  childs_edge[0][0] = 0; childs_edge[0][1] = 1;
  childs_edge[1][0] = 1; childs_edge[1][1] = 2;
  childs_edge[2][0] = 2; childs_edge[2][1] = 3;
  childs_edge[3][0] = 3; childs_edge[3][1] = 0;
  childs_edge[4][0] = 4; childs_edge[4][1] = 5;
  childs_edge[5][0] = 5; childs_edge[5][1] = 6;
  childs_edge[6][0] = 6; childs_edge[6][1] = 7;
  childs_edge[7][0] = 7; childs_edge[7][1] = 4;
  childs_edge[8][0] = 0; childs_edge[8][1] = 4;
  childs_edge[9][0] = 1; childs_edge[9][1] = 5;
  childs_edge[10][0] = 2; childs_edge[10][1] = 6;
  childs_edge[11][0] = 3; childs_edge[11][1] = 7;

  //  vicf[face][child] = local number of face vertex
  vicf[0][0] = 2; vicf[0][1] = 3; vicf[0][2] = 0; vicf[0][3] = 1;
  vicf[1][0] = 6; vicf[1][1] = 2; vicf[1][2] = 1; vicf[1][3] = 5;
  vicf[2][0] = 6; vicf[2][1] = 7; vicf[2][2] = 3; vicf[2][3] = 2;
  vicf[3][0] = 3; vicf[3][1] = 7; vicf[3][2] = 4; vicf[3][3] = 0;
  vicf[4][0] = 1; vicf[4][1] = 0; vicf[4][2] = 4; vicf[4][3] = 5;
  vicf[5][0] = 7; vicf[5][1] = 6; vicf[5][2] = 5; vicf[5][3] = 4;

  //  vice[edge][child] = local number of edge vertex
  vice[0][0] = 1; vice[0][1] = 0; 
  vice[1][0] = 2; vice[1][1] = 1; 
  vice[2][0] = 3; vice[2][1] = 2; 
  vice[3][0] = 0; vice[3][1] = 3; 
  vice[4][0] = 5; vice[4][1] = 4; 
  vice[5][0] = 6; vice[5][1] = 5; 
  vice[6][0] = 7; vice[6][1] = 6; 
  vice[7][0] = 4; vice[7][1] = 7; 
  vice[8][0] = 4; vice[8][1] = 0; 
  vice[9][0] = 5; vice[9][1] = 1; 
  vice[10][0] = 6; vice[10][1] = 2; 
  vice[11][0] = 7; vice[11][1] = 3; 

  // gibt die bzgl. der Kinder lokale vertex number des innersten Knoten
  cell_midpoint[0] = 6;
  cell_midpoint[1] = 7;
  cell_midpoint[2] = 4;
  cell_midpoint[3] = 5;
  cell_midpoint[4] = 2;
  cell_midpoint[5] = 3;
  cell_midpoint[6] = 0;
  cell_midpoint[7] = 1;

  // gibt die lokalen inneren face nummern der Kinder
  ieoc[0][0] = 1; ieoc[0][1] = 5; ieoc[0][2] = 2;
  ieoc[1][0] = 3; ieoc[1][1] = 5; ieoc[1][2] = 2;
  ieoc[2][0] = 3; ieoc[2][1] = 4; ieoc[2][2] = 5;
  ieoc[3][0] = 1; ieoc[3][1] = 4; ieoc[3][2] = 5;

  ieoc[4][0] = 0; ieoc[4][1] = 1; ieoc[4][2] = 2;
  ieoc[5][0] = 0; ieoc[5][1] = 3; ieoc[5][2] = 2;
  ieoc[6][0] = 3; ieoc[6][1] = 4; ieoc[6][2] = 0;
  ieoc[7][0] = 1; ieoc[7][1] = 4; ieoc[7][2] = 0;

  // gibt die lokalen nummern der Kinder zu einer inneren face 
  coif[0][0] = 0; coif[0][1] = 1;
  coif[1][0] = 1; coif[1][1] = 2;
  coif[2][0] = 2; coif[2][1] = 3;
  coif[3][0] = 3; coif[3][1] = 0;
  coif[4][0] = 4; coif[4][1] = 5;
  coif[5][0] = 5; coif[5][1] = 6;
  coif[6][0] = 6; coif[6][1] = 7;
  coif[7][0] = 7; coif[7][1] = 4;
  coif[8][0] = 0; coif[8][1] = 4;
  coif[9][0] = 1; coif[9][1] = 5;
  coif[10][0] = 2; coif[10][1] = 6;
  coif[11][0] = 3; coif[11][1] = 7;

  fixarray<6,int> oppositeface;
  
  oppositeface[0] = 5;
  oppositeface[1] = 3;
  oppositeface[2] = 4;
  oppositeface[3] = 1;
  oppositeface[4] = 2;
  oppositeface[5] = 0;

  // lcfif[e][i] gibt zur e-ten inneren face die lokale face
  //             nummer des i-ten hex dass daran liegt an

  lcfif[0][0] = 1;
  lcfif[1][0] = 2;
  lcfif[2][0] = 3; 
  lcfif[3][0] = 4; 
  lcfif[4][0] = 1; 
  lcfif[5][0] = 2; 
  lcfif[6][0] = 3; 
  lcfif[7][0] = 4; 
  lcfif[8][0] = 5; 
  lcfif[9][0] = 5; 
  lcfif[10][0] = 5;
  lcfif[11][0] = 5;

  for (int i=0; i<12; i++)
    {
      lcfif[i][1] = oppositeface[ lcfif[i][0] ];
    }

  edgeofface[0][0] = 0;
  edgeofface[0][1] = 1;
  edgeofface[0][2] = 2;
  edgeofface[0][3] = 3;

  edgeofface[1][0] = 9;
  edgeofface[1][1] = 5;
  edgeofface[1][2] = 10;
  edgeofface[1][3] = 1;

  edgeofface[2][0] = 2;
  edgeofface[2][1] = 10;
  edgeofface[2][2] = 6;
  edgeofface[2][3] = 11;

  edgeofface[2][0] = 10;  // Reihenfolge geaendert 09.Jan.2002, Malte
  edgeofface[2][1] = 6;
  edgeofface[2][2] = 11;
  edgeofface[2][3] = 2;

  edgeofface[3][0] = 11;  // Reihenfolge geaendert !!!
  edgeofface[3][1] = 7;
  edgeofface[3][2] = 8;
  edgeofface[3][3] = 3;

  edgeofface[4][0] = 8;  // Reihenfolge geaendert !!!
  edgeofface[4][1] = 4;
  edgeofface[4][2] = 9;
  edgeofface[4][3] = 0;

  edgeofface[5][0] = 7;  // Reihenfolge geaendert !!!
  edgeofface[5][1] = 6;
  edgeofface[5][2] = 5;
  edgeofface[5][3] = 4;

  // for Q2 hanging nodes on patched faces

  hnpf[0][0] = 0; hnpf[1][0] = 2; hnpf[2][0] = 6; hnpf[3][0] = 8;
  hnpf[0][1] = 1; hnpf[1][1] = 5; hnpf[2][1] = 3; hnpf[3][1] = 7;
  hnpf[0][2] = 2; hnpf[1][2] = 8; hnpf[2][2] = 0; hnpf[3][2] = 6;
  hnpf[0][3] = 3; hnpf[1][3] = 1; hnpf[2][3] = 7; hnpf[3][3] = 5;
  hnpf[0][4] = 4; hnpf[1][4] = 4; hnpf[2][4] = 4; hnpf[3][4] = 4;
  hnpf[0][5] = 5; hnpf[1][5] = 7; hnpf[2][5] = 1; hnpf[3][5] = 3;
  hnpf[0][6] = 6; hnpf[1][6] = 0; hnpf[2][6] = 8; hnpf[3][6] = 2;
  hnpf[0][7] = 7; hnpf[1][7] = 3; hnpf[2][7] = 5; hnpf[3][7] = 1;
  hnpf[0][8] = 8; hnpf[1][8] = 6; hnpf[2][8] = 2; hnpf[3][8] = 0;
}

/*----------------------------------------------------------------------*/
/* laedt in E zum Hex q und einer face F die e-te Kante                 */

int HexLawAndOrder::EdgeVertexOfFace(const Hex& q, const FaceVector& F, int e) const
{
  int face = local_face(q,F);
  int le = edgeofface[face][e];
  return edge_vertex(q,le);
}

/*----------------------------------------------------------------------*/

int HexLawAndOrder::LoadEdgeOfFace(const Hex& q, const FaceVector& F, int e, EdgeVector& E
) const
{
  int face = local_face(q,F);
  int edge = edgeofface[face][e];
  globalvertices_of_edge(q,E,edge);
  return edge;
}

/*----------------------------------------------------------------------*/

// gibt ne globale face Nummer zurueck
int HexLawAndOrder::InnerEdge(const Hex& h, int i) const
{
  int j = lcfif[i][0];
  int k = coif[i][0];

  return hexs[ h.child(k) ].edge(j);
}

/*----------------------------------------------------------------------*/

int HexLawAndOrder::GlobalInnerFace(int c, int i) const
{
  const Hex& h = hexs[c];
  int ie = InnerEdge(h,i);
  return h.edge(ie);
}

/*----------------------------------------------------------------------*/

// gibt ne globale face Nummer zurueck
int HexLawAndOrder::GlobalChildFace(const FaceVector& face, int q, int j) const
{
  int ledge = local_face_index(q,face);
  int ic    = ChildsOfFace(ledge,j);
  int child = hexs[q].child(ic);
  int iedge = ChildFace(ledge);
  
  return hexs[child].edge(iedge);
}

/*---------------------------------------------------*/

// laedt vertex indices in facevector
void HexLawAndOrder::GetFace(FaceVector& face, int h, int e) const
{
  const Hex& H = hexs[h];
  for (int i=0; i<4; i++)
    {
      face[i] = H[ lvf[e][i] ];
    }
}

/*---------------------------------------------------*/

// laedt vertex indices in edgevector
void HexLawAndOrder::globalvertices_of_edge(const Hex& q, fixarray<2,int>& f, int edge) const
{
  local_edge_index(f, edge);

  f[0] = q.vertex(f[0]);
  f[1] = q.vertex(f[1]);
}

/*---------------------------------------------------*/

int HexLawAndOrder::GetVertexOfEdge(int iq, const fixarray<2,int>& edge) const
{
  assert(iq>=0);

  const Hex& H = hexs[iq];
  
  fixarray<2,int> v;
  EdgeArray<2> Edge(edge);

  for (int e=0; e<12; e++)
    {
      globalvertices_of_edge(H,v,e);
      if (Edge==v) 
	{
	  int c0 = childs_edge[e][0];
	  int ge = hexs[H.child(c0)].vertex(vice[e][0]);
	  
	  return ge;
	}
    }
  abort();
}

/*---------------------------------------------------*/

void HexLawAndOrder::globalvertices_of_face(const Hex& q, fixarray<4,int>& f, int face) const
{
  local_face_index(f, face);

  f[0] = q.vertex(f[0]);
  f[1] = q.vertex(f[1]);
  f[2] = q.vertex(f[2]);
  f[3] = q.vertex(f[3]);
}

/*---------------------------------------------------*/

pair<int,int> HexLawAndOrder::GetChildFaces(const FaceVector& bigedge, 
					    int bighex, int i) const
{
  // bigedge sind die vier vertex no der grossen face
  // bighex ist die no des vater hex
  // i ist die no der kleinen Kinderface (0<=i<=4)

  int bigeind = local_face_index(bighex,bigedge);

  assert(bigeind>=0);
  assert(bigeind< 6);

//   if ((bigeind>=6) || (bigeind<0))
//     {
//       cout << "\nHexLawAndOrder::GetChildFaces bigeind ";
//       cout << bighex << " [" << bigedge << "] -> " << bigeind << endl; 
//       cout << "Hex was " << hexs[bighex] << endl;
//       return make_pair(-1,-1);
//       //exit(1);
//     }

  int iedge   = ChildFace(bigeind);
  int ic      = ChildsOfFace(bigeind,i);
  int q       = hexs[bighex].child(ic);
  
  return make_pair(q,iedge);
}

/*----------------------------------------------------------------------*/

int HexLawAndOrder::local_face_index(int q, const FaceVector& face) const
{
  EdgeArray<4> f(face);
  EdgeArray<4> g(face);
  for (int i=0; i<6; i++)
    {
      global_face_unsorted(g,hexs[q],i);
      if (f==g) return i;
    }
  return -1;
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::local_edge_index(EdgeVector& v, int edge) const
{
  v = lve[edge];
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::local_face_index(FaceVector& v, int face) const
{
  v = lvf[face];
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::global_face_unsorted(FaceVector& faceglob, const Hex& h,
					  int face) const
{
  FaceVector faceloc;
  local_face_index(faceloc,face);
  h.vertex_loc2glob(faceglob,faceloc);
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::global_edge_unsorted(EdgeVector& edgeglob, const Hex& h,
					  int edge) const
{
  EdgeVector edgeloc;
  local_edge_index(edgeloc,edge);
  h.vertex_loc2glob(edgeglob,edgeloc);
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::fill_corner_vertex_in_childs(const Hex& f) const
{
  hexs[f.child(0)].vertex(0) = f.vertex(0);
  hexs[f.child(1)].vertex(1) = f.vertex(1);
  hexs[f.child(2)].vertex(2) = f.vertex(2);
  hexs[f.child(3)].vertex(3) = f.vertex(3);
  hexs[f.child(4)].vertex(4) = f.vertex(4);
  hexs[f.child(5)].vertex(5) = f.vertex(5);
  hexs[f.child(6)].vertex(6) = f.vertex(6);
  hexs[f.child(7)].vertex(7) = f.vertex(7);
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::fill_middle_vertex_in_childs(const Hex& f,
						  int number) const  
{
  hexs[f.child(0)].vertex(6) = number;
  hexs[f.child(1)].vertex(7) = number;
  hexs[f.child(2)].vertex(4) = number;
  hexs[f.child(3)].vertex(5) = number;
  hexs[f.child(4)].vertex(2) = number;
  hexs[f.child(5)].vertex(3) = number;
  hexs[f.child(6)].vertex(0) = number;
  hexs[f.child(7)].vertex(1) = number;
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::fill_face_vertex_in_childs(const Hex& f, 
						int face, int number) const
{
  int c0 = childs_face[face][0];
  int c1 = childs_face[face][1];
  int c2 = childs_face[face][2];
  int c3 = childs_face[face][3];

  hexs[f.child(c0)].vertex(vicf[face][0]) = number;
  hexs[f.child(c1)].vertex(vicf[face][1]) = number;
  hexs[f.child(c2)].vertex(vicf[face][2]) = number;
  hexs[f.child(c3)].vertex(vicf[face][3]) = number;
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::fill_edge_vertex_in_childs(const Hex& f, 
						int edge, int number) const
{
  int c0 = childs_edge[edge][0];
  int c1 = childs_edge[edge][1];

  hexs[f.child(c0)].vertex(vice[edge][0]) = number;
  hexs[f.child(c1)].vertex(vice[edge][1]) = number;
}

/*-------- for regular one ---------------------------------------------*/

void HexLawAndOrder::childs_of_global_face(FaceVector& child, const Hex& f, 
					   const FaceVector& globalface) const
{
  int localface = local_face(f,globalface);
  childs_of_face(child,f,localface);
}

/*---------------------------------------------------------------------*/

int HexLawAndOrder::TestFaceOfOneChild(const Hex& H, const FaceVector& F) const
{
  int lf = local_face(H,F);

  int childhex = 0;
  int child = childs_face[lf][childhex];
  int gchild = H.child(child);
  int lcf = ChildFace(lf);
  int res = hexs[gchild].edge(lcf);
  return res;
}

/*---------------------------------------------------------------------*/

int HexLawAndOrder::local_face(const Hex& f, 
			       const FaceVector& globalface) const
{
  // gives the local face number of a face

  EdgeArray<4> li(globalface);  // wird gleich wieder ueberschrieben !!
  for (int i=0; i<4; i++)
    {
      li[i]= f.global2local(globalface[i]);
    }
  for (int i=0; i<6; i++)
    {
      if (li==lvf[i]) return i;
    }
  abort();
}

/*---------------------------------------------------------------------*/

void HexLawAndOrder::childs_of_face(FaceVector& child, const Hex& f, 
				    int face) const
{
  child[0] = f.child(childs_face[face][0]);
  child[1] = f.child(childs_face[face][1]);
  child[2] = f.child(childs_face[face][2]);
  child[3] = f.child(childs_face[face][3]);
}

/*----------------------------------------------------------------------*/

int HexLawAndOrder::face_vertex(const Hex& f, int face) const
{
  if (!f.sleep()) return -1;

  int lic = childs_face[face][0];
  int liv = vicf[face][0];
  int gic = f.child(lic);

  return hexs[gic].vertex(liv);
}

/*----------------------------------------------------------------------*/

// fuer boundary newton
void HexLawAndOrder::LoadEdgeVerticesOfFace(const Hex& f, int face, FaceVector& dst) const
{
  for (int i=0; i<4; i++) 
    {
      int edge = edgeofface[face][i];
      dst[i] = edge_vertex(f,edge);
    }
}

/*----------------------------------------------------------------------*/

// fuer boundary newton
void HexLawAndOrder::LoadFaceVertices(const Hex& f, fixarray<6,int>& dst) const
{
  for (int i=0; i<6; i++) 
    {
      dst[i] = face_vertex(f,i);
    }
}

/*----------------------------------------------------------------------*/

int HexLawAndOrder::edge_vertex(const Hex& f, int edge) const
{
  int lic = childs_edge[edge][0]; // local edgenumber of corresponding child
  int liv = vice[edge][0];
  int gic = f.child(lic);

  return hexs[gic].vertex(liv);
}

/*----------------------------------------------------------------------*/

int HexLawAndOrder::middle_vertex(const Hex& f) const
{
  return hexs[f.child(0)].vertex(cell_midpoint[0]);
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::load_face(FaceVector& v, const Hex& f, int face) const
{
  for (int i=0; i<4; i++)
    {
      v[i] = f.vertex(lvf[face][i]);
    }
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::GetGlobalOuterFaceOfChild
(FaceVector& face, const Hex& f, int c, int e) const
{
  int childface = 0;
  const Hex& child = hexs[f.child(c)];
  load_face(face,child,childface);
}

/*----------------------------------------------------------------------*/

void HexLawAndOrder::globalfacechildren_of_father(vector<FaceVector>& faces,
						   const Hex& f) const
{
  size_t n = 24;
  faces.resize(n);

  int k = 0;
  for (int i=0; i<8; i++)
    {
      for (int e=0; e<3; e++)
	{
	  GetGlobalOuterFaceOfChild(faces[k++],f,i,e);
	}
    }      
}

/*---------------------------------------------------*/

fixarray<9,int> HexLawAndOrder::PatchVerticesOfFace(int h, int face) const
{
  fixarray<4,int> f;
  GetFace(f,h,face);
  //  globalvertices_of_face(q,f,face);

  fixarray<9,int> F;
  F[0] = f[3];
  F[2] = f[2];
  F[6] = f[0];
  F[8] = f[1];

  const Hex& q = hexs[h];
  F[4] = face_vertex(q,face);

  F[1] = EdgeVertexOfFace(q,f,2);
  F[3] = EdgeVertexOfFace(q,f,3);
  F[5] = EdgeVertexOfFace(q,f,1);
  F[7] = EdgeVertexOfFace(q,f,0);

  return F;
}

/*---------------------------------------------------*/

fixarray<9,int> HexLawAndOrder::GiveOrdering(const fixarray<9,int>& F, const Hex& h) const
{
  int found = -1;
  fixarray<4,int> kk; kk[0] = 0; kk[1] = 2; kk[2] = 6; kk[3] = 8;
  for (int j=0; (j<4) && (found<0); j++)
    {
      int fk = F[kk[j]];
      for (int i=0; (i<8) && (found<0); i++)
	{
	  if (h[i]==fk)  found = j;
	}
    }
  assert(found>=0); 
  return hnpf[found];
}
}
