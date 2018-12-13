/**
*
* Copyright (C) 2004, 2008, 2011 by the Gascoigne 3D authors
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


#include  "facemanager.h"
#include  "vecalgo.h"
#include  "hangfacesort.h"
#include  "giota.h"

#include  <map>
#ifdef __OLDCOMPILER__
#include  <hash_map>
#define HANGMAP  hash_map<EdgeArray<4>,int,EdgeHash>
#else
#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HANGMAP   std::tr1::unordered_map<EdgeArray<4>,int,EdgeHash> 
#else
#include  <ext/hash_map>
#define HANGMAP  __gnu_cxx::hash_map<EdgeArray<4>,int,EdgeHash> 
#endif
#endif

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne
{
FaceManager::FaceManager(vector<Edge>& e, vector<Hex>& q, 
			 const IntVector& con, IntVector& eon) :
  edges(e), hexs(q), co2n(con), eo2n(eon), HexLaO(q) {}

/*---------------------------------------------------*/

void FaceManager::InitFaces()
{
  HANGMAP H;

  FaceVector e;

  for (int i=0; i<hexs.size(); i++)
    {
      for (int j=0; j<6; j++)
	{
	  HexLaO.global_face_unsorted(e,hex(i),j);

	  HANGMAP::iterator yes =  H.find(e);

	  if (yes!=H.end())
	    {
	      int k = yes->second;
	      edges[k].slave() = i;
	      edges[k].LocalSlaveIndex()=j;
	      hexs[i].edge(j) = k;
	      H.erase(yes);
	    }
	  else
	    {
	      Edge E(i,j);
	      int n = edges.size();
	      edges.push_back(E);
	      H.insert(make_pair(e,n));
	      hexs[i].edge(j) = n;
	    }
	}
    }
  SortHangings();
}

/*---------------------------------------------------*/

void FaceManager::DeleteFaces()
{
  compress(edges,eo2n);

  for (int i=0; i<hexs.size(); i++)
    {
      Hex& H = hexs[i];
      if (co2n[i]<0)
	{
	  H.edges() = -1;
	}
      else
	{
	  for (int e=0; e<6; e++)
	    {
	      int edge = H.edge(e);
	      int newe = eo2n[ edge ];

	      assert(newe>=0);

	      H.edge(e) = newe;
	    }
	}
    }
}

/*---------------------------------------------------*/

void FaceManager::Update()
{
  for (int i=0; i<edges.size(); i++)
    {
      int m = edges[i].master();
      int s = edges[i].slave();

      assert(m>=0);

      int nm = co2n[m];
      if (nm>=0)
	{
	  edges[i].master() = nm;
	  if (s>=0)
	    {
	      edges[i].slave() = co2n[s];
	    }
	}
      else
	{
	  if (s<0)
	    {
	      edges[i].master() = -1;
	      SwappedEdge.push_back(i);
	      continue;
	    }
	  edges[i].swapping(co2n[s]);
	}
    }
}

/*---------------------------------------------------*/

void FaceManager::Build(const IntSet& CellRefList,
			HangContainer3d& hangset)
{
  SwappedEdge.resize(0);

  Update     ();
  InnerFaces (CellRefList);
  OuterFaces (hangset);
  OldHangings(hangset,CellRefList);
  SwappedFaces(); 
  NeighbourTester();
  SortHangings();
  Check(hangset);
}

/*---------------------------------------------------*/

void FaceManager::InnerFaces(const IntSet& CellRefList)
{
  int n   = edges.size();
  int nv1 = CellRefList.size();

  edges.reserve(n+12*nv1);
  edges.resize(n+12*nv1);

  IntSet::const_iterator cp;

  for (cp=CellRefList.begin(); cp!=CellRefList.end(); cp++)
    {
      int f = co2n[*cp];

      for (int e=0; e<12; e++)
	{
	  int icl = HexLaO.ChildOfInnerFace(e,0);
	  int ic  = hexs[f].child(icl);
	  int ie  = HexLaO.LocalChildFaceOfInnerFace(e,0); 
	  
	  hexs[ic].edge(ie) = n;

	  Edge E(ic,ie);

	  int icl2 = HexLaO.ChildOfInnerFace(e,1);
	  int ic2  = hexs[f].child(icl2);
	  int ie2  = HexLaO.LocalChildFaceOfInnerFace(e,1);
	  hexs[ic2].edge(ie2) = n;

	  E.slave()           = ic2;
	  E.LocalSlaveIndex() = ie2;

	  edges[n] = E;
	  n++;
	}
    }
}

/*---------------------------------------------------*/

void FaceManager::Check(const HangContainer3d& hangset) const
{
  HangList<4>::const_iterator p = hangset.FaceNotMore().begin();

  int ok = 1;
  for (;p!=hangset.FaceNotMore().end(); p++)
    {
      int hang = p->second.hanging();
      int M    = p->second.rneighbour();
      int S    = p->second.cneighbour();

      if (hang<0) continue;
      if (M<0)    continue;
      if (S<0)    continue;   

      const Hex& HM = hexs[M];
      const Hex& HS = hexs[S];

      if (!HM.sleep())
	{
	  cout << "Hex master has no childs " << M << endl;
	  cout << "face " << p->first << endl;
	  ok = 0;
	}
      if (!HS.sleep())
	{
	  cout << "Hex slave has no childs " << S << endl;
	  cout <<  "face " << p->first << endl;
	  cout << " slave " << HS << endl;
	  ok = 0;
	}
    }
  assert(ok);
}

/*---------------------------------------------------*/

void FaceManager::SortHangings()
{
  // edges with hanging nodes swapped to the end of list

  vector<int> perm(edges.size());
  iota(perm.begin(),perm.end(),0);
  stable_sort(perm.begin(),perm.end(),HangFaceSort(*this));
  stable_sort(edges.begin(),edges.end(),HangFaceSort2(*this));
 
  int i=edges.size();
  while(EdgeIsHanging(i-1)) i--;

  vector<int> permi(perm.size());
  for(int i=0;i<perm.size();i++) permi[perm[i]] = i;

  for(int i=0;i<hexs.size();i++)
    {
      for(int ii=0;ii<6;ii++)  
	{
	  hexs[i].edge(ii) = permi[hex(i).edge(ii)];
	}
    }
  for (int i=0; i<eo2n.size(); i++)
    {
      if (eo2n[i]>=0) 
	{
	  eo2n[i] = permi[eo2n[i]]; 
	}
    }

  // master of each edge is allways the coarser hex

  for(int i=0;i<edges.size();i++)
    {
      int m = edges[i].master();
      int s = edges[i].slave();
      if(s<0) continue;
      if (hex(m).sleep() && !hex(s).sleep())
	{ 
	  //cout << "Swapping " << m << " " << s << endl;
	  swap(edges[i].master(),edges[i].slave());
	  swap(edges[i].LocalMasterIndex(),edges[i].LocalSlaveIndex());
	}
    }
}

/*---------------------------------------------------*/

bool FaceManager::EdgeIsHanging(int e) const
{
  int m = edges[e].master();
  int s = edges[e].slave();
  if (s<0) return 0;
  if (hex(m).sleep() && !hex(s).sleep() ) return 1;
  if (hex(s).sleep() && !hex(m).sleep() ) return 1;
  return 0;
}

/*---------------------------------------------------*/

bool FaceManager::EdgeIsHanging(const Edge& e) const
{
  int m = e.master();
  int s = e.slave();
  if (s<0) return 0;
  if (hex(m).sleep() && !hex(s).sleep() ) return 1;
  if (hex(s).sleep() && !hex(m).sleep() ) return 1;
  return 0;
}

/*---------------------------------------------------*/

void FaceManager::LoadFaceElimination(IntVector& edel, 
				      const IntSet& CellCoarseList,
				      const HangContainer3d& hangset) const
{
  edel.resize(12*CellCoarseList.size());

  int n = 0;
  
  IntSet::const_iterator cp;
  
  for (cp=CellCoarseList.begin(); cp!=CellCoarseList.end(); cp++)
    {
      if (*cp<0) continue;

      const Hex& H = hexs[*cp];
      for (int f=0; f<12; f++)
 	{
	  edel[n++] = HexLaO.InnerEdge(H,f);
	}
    }
  HangList<4>::const_iterator p = hangset.FaceDeleting().begin();

  edel.resize(n+4*hangset.nFaceVertexesToBeDeleted());

  for (;p!=hangset.FaceDeleting().end(); p++)
    {
      const FaceVector& F = p->first;
      int h = p->second.rneighbour();
      for (int i=0; i<4; i++)
 	{
	  edel[n++] = HexLaO.GlobalChildFace(F,h,i);
  	}
    }
  edel.resize(n);

  sort(edel.begin(),edel.end());
}

/*---------------------------------------------------*/

void FaceManager::NeighbourTester() const
{
  IntVector x(edges.size());
  for (int i=0; i<hexs.size(); i++)
    {
      for (int e=0; e<6; e++)
	{
	  int edge = hexs[i].edge(e);
	  x[edge]++;
	}
    }
  for (int i=0; i<x.size(); i++)
    {
      if (x[i]>2)
	{
	  cout << "BSE Test " << i << " " << x[i] << endl;
	}
    }
  for (int i=0; i<hexs.size(); i++)
    {
      for (int e=0; e<6; e++)
	{
	  int edge = hexs[i].edge(e);

	  if (edge<0) continue;
	  const Edge& E = edges[edge];
	  int m = E.master();
	  int s = E.slave();
	  int nachbar = -10;
	  if      (m==i) { nachbar = s; }
	  else if (s==i) { nachbar = m; }

	  if (nachbar==-10) 
	    {
	      cout << "BSE Test " << i << " " << e << " edgenumber " << edge << endl;
	    }
	}
    }
}

/*---------------------------------------------------*/

void FaceManager::FillNeighbourFaces(const Hex& HM, const Hex& HS,
				     const FaceVector& Face)
{
  int MF = HexLaO.local_face(HM,Face);
  int SF = HexLaO.local_face(HS,Face);

  FaceVector m,s;
  HexLaO.childs_of_face(m,HM,MF);
  HexLaO.childs_of_face(s,HS,SF);

  int mf = HexLaO.ChildFace(MF);
  int sf = HexLaO.ChildFace(SF);

  //  cout <<  endl << "FACE " << Face << endl;
  for (int i=0; i<4; i++)
    {
      int master = m[i];

      FaceVector face;
      HexLaO.GetFace(face,master,mf);

      //      cout << master <<" master child face " << face << endl;

      int slave = -1; 
      int sj    = -1;
      for (int j=0; (j<4)&&(sj<0); j++)
	{
	  slave = s[j];
	  sj = HexLaO.local_face_index(slave,face);
	}
      assert(sj>=0);
      assert(sj==sf);

            Hex& hs = hexs[ slave ];
      const Hex& hm = hexs[ master];
      int e = hm.edge(mf);

      assert(hs.edge(sj)==-1);

      hs.edge(sj) = e;
      Edge& E = edges[e];

      if (E.master()!=master)
	{
	  assert(E.slave()==master);

	  E.master() = master;
	  E.LocalMasterIndex() = E.LocalSlaveIndex();	  
	}
      E.slave() = slave;
      E.LocalSlaveIndex() = sf;
    }
}

/*---------------------------------------------------*/

void FaceManager::OuterFaces(const HangContainer3d& hangset)
{
  int n = edges.size();

  edges.reserve(n+4*hangset.FaceCreating().size());
  edges.resize (n+4*hangset.FaceCreating().size());

  HangList<4>::const_iterator p = hangset.FaceCreating().begin();

  // neue Kinder-faces erzeugen
  for (;p!=hangset.FaceCreating().end(); p++)
    {
      int rneigh = p->second.rneighbour();

      for (int i=0; i<4; i++)
	{
	  pair<int,int> cp = HexLaO.GetChildFaces(p->first,rneigh,i);

	  int cellindex = cp.first;
	  int edgeindex = cp.second;

	  hexs[cellindex].edge(edgeindex) = n;

	  edges[n++] = Edge(cellindex,edgeindex);
	}
    }
  p = hangset.FaceCreating().begin();

  for (;p!=hangset.FaceCreating().end(); p++)
    {
      int M = p->second.rneighbour();
      int S = p->second.cneighbour();
      // wenn neue Kinder-faces nicht hangen, dann Nachbarn eintragen
      if (S>=0)
	{
	  const Hex& HM = hexs[M];
	  const Hex& HS = hexs[S];
	  if (HS.sleep())
	    {
	      FillNeighbourFaces(HM,HS,p->first);
	    }
	}
    }
}

/*---------------------------------------------------*/

void FaceManager::OldHangings(HangContainer3d& hangset,
			      const IntSet& CellRefList)
{
  HangList<4>::iterator p = hangset.FaceNotMore().begin();

  for (;p!=hangset.FaceNotMore().end(); p++)
    {
      int hang = p->second.hanging();
      int M    = p->second.rneighbour();
      int S    = p->second.cneighbour();

      if (hang<0) continue;
      if (M<0)    continue;
      if (S<0)    continue;   

      const Hex& HM = hexs[M];
      const Hex& HS = hexs[S];

      const FaceVector& F = p->first;

//       cout << M << " Master " << HM;
//       for (int i=0; i<HS.nchilds(); i++)
// 	{
// 	  int j = HM.child(i);
// 	  cout << j << "child " << hexs[j];
// 	}
//       cout << S << " Slave " << HS << endl;
	  
//       for (int i=0; i<HS.nchilds(); i++)
// 	{
// 	  int j = HS.child(i);
// 	  cout << j << " hild " << hexs[j];
// 	}
      
      int em = HexLaO.TestFaceOfOneChild(HM,F);
      int es = HexLaO.TestFaceOfOneChild(HS,F);

      if (es<0)
	{
	  FillNeighbourFaces(HM,HS,F);
	}
      else if (em<0)
	{
	  FillNeighbourFaces(HS,HM,F);
	}
    }
}

/*---------------------------------------------------*/

void FaceManager::SwappedFaces()
{
  int n = 0;
  int m = 0;
  for (int i=0; i<hexs.size(); i++)
    {
      const Hex& q = hexs[i];
      for (int e=0; e<6; e++)
	{
	  if (q.edge(e)<0) m++;
	}
    }
  assert(m==SwappedEdge.size());

  for (int i=0; i<hexs.size(); i++)
    {
      Hex& q = hexs[i];
      if (q.sleep()) continue;
      for (int e=0; e<6; e++)
	{
	  if (q.edge(e)<0)
	    {
	      int ei = SwappedEdge[n++];
	      q.edge(e) = ei;

	      edges[ei].setmaster(i,e);
	    }
	}
    }
}
}

#undef HANGMAP
