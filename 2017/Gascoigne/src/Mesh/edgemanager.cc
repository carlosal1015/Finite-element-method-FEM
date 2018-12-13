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


#include  "edgemanager.h"
#include  "vecalgo.h"
#include  "hangsort.h"
#include  "giota.h"

#include  <map>
#ifdef __OLDCOMPILER__
#include  <hash_map>
#define HANGMAP  hash_map<EdgeArray<2>,int,EdgeHash>
#else
#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HANGMAP std::tr1::unordered_map<EdgeArray<2>,int,EdgeHash> 
#else
#include  <ext/hash_map>
#define HANGMAP __gnu_cxx::hash_map<EdgeArray<2>,int,EdgeHash> 
#endif
#endif

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne
{
EdgeManager::EdgeManager(vector<Edge>& e, vector<Quad>& q, 
			 const IntVector& con, IntVector& eon) :
  edges(e), quads(q), co2n(con), eo2n(eon), QuadLaO(q) {}

/*---------------------------------------------------*/

void EdgeManager::BSETest() const
{
  cout << "BSE Tester:\n";
  IntVector x(edges.size());
  for (int i=0; i<quads.size(); i++)
    {
      for (int e=0; e<4; e++)
	{
	  int edge = quads[i].edge(e);
	  x[edge]++;
	}
    }
  for (int i=0; i<x.size(); i++)
    {
      if (x[i]>2)
	{
	  cout << "problem 1 in edge " << i << " " << x[i] << endl;
	}
    }
  for (int i=0; i<quads.size(); i++)
    {
      for (int e=0; e<4; e++)
	{
	  int edge = quads[i].edge(e);

	  if (edge<0) continue;
	  const Edge& E = edges[edge];
	  int m = E.master();
	  int s = E.slave();
	  int ml = E.LocalMasterIndex();
	  int sl = E.LocalSlaveIndex();
	  int nachbar = -10;
	  int f = -1;
	  if      (m==i) { nachbar = s; f = sl;}
	  else if (s==i) { nachbar = m; f = ml;}

	  if (nachbar==-10) 
	    {
	      cout << "problem 2 in edge " << i << " " << e << " edgenumber " << edge << " f= " << f << endl;
	    }
	}
    }
}

/*---------------------------------------------------*/

void EdgeManager::LoadEdgeElimination(IntVector& edel, 
				      const IntSet& CellCoarseList,
				      const HangContainer2d& hangset) const
{
  edel.resize(4*CellCoarseList.size()+2*hangset.NToBeDeleted());

  int n = 0;
  
  IntSet::const_iterator cp;

  for (cp=CellCoarseList.begin(); cp!=CellCoarseList.end(); cp++)
    {
      if (*cp<0) continue;

      for (int e=0; e<4; e++)
	{
	  edel[n++] = QuadLaO.GlobalInnerEdge(*cp,e);
	}
    }
  HangList<2>::const_iterator p = hangset.Deleting().begin();

  for (;p!=hangset.Deleting().end(); p++)
    {
      for (int i=0; i<2; i++)
	{
	  edel[n++] = QuadLaO.GlobalChildEdge(p->first,p->second.rneighbour(),i);
	}
    }
  edel.resize(n);
}

/*---------------------------------------------------*/

void EdgeManager::Build(const IntSet& CellRefList,
		       HangContainer2d& hangset)
{
  IntVector SwappedEdge;

  Update     ();
  InnerEdges (CellRefList);
  //BSETest();
  OuterEdges (hangset);
  OldHangings(hangset,CellRefList);
  SwappedEdges();
  NeighbourTester();
  SortHangings();
}

/*---------------------------------------------------*/

void EdgeManager::Update()
{
  for (int i=0; i<edges.size(); i++)
    {
      int m = edges[i].master();
      int s = edges[i].slave();
      if (m<0)
	{
	  cerr << "EdgeManager::update()" << endl;
	  cerr << "Master negativ " << i << " " << m << endl;
	  abort();
	}
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

void EdgeManager::InnerEdges(const IntSet& CellRefList)
{
  int n   = edges.size();
  int nv1 = CellRefList.size();
  
  edges.reserve(n+4*nv1);
  edges.resize(n+4*nv1);

  IntSet::const_iterator cp;

  for (cp=CellRefList.begin(); cp!=CellRefList.end(); cp++)
    {
      int f = co2n[*cp];

      for (int e=0; e<4; e++)
	{
	  int ic = quad(f).child(e);
	  int ie = QuadLaO.InnerEdgeOfChild(e,0);
	  quad(ic).edge(ie) = n;

	  int oe = QuadLaO.OuterEdgeOfChild(e,0);
	  quads[ic].edge(oe) = -1;
	  oe = QuadLaO.OuterEdgeOfChild(e,1);
	  quad(ic).edge(oe) = -1;

	  if (ic<0)
	    {
	      cout << "master = " << ic << endl;
	    }

	  Edge E(ic,ie);

	  ic = quad(f).child((e+1)%4);
	  ie = QuadLaO.InnerEdgeOfChild((e+1)%4,1);

	  quad(ic).edge(ie) = n;

	  E.slave()           = ic;
	  E.LocalSlaveIndex() = ie;

	  edges[n] = E;
	  n++;
	}
    }
}

/*---------------------------------------------------*/

void EdgeManager::OuterEdges(const HangContainer2d& hangset)
{
  int n = edges.size();
  int nn = n+2*hangset.NToBeCreated();

  edges.reserve(nn);
  edges.resize (nn);

  HangList<2>::const_iterator p = hangset.Creating().begin();

  for (;p!=hangset.Creating().end(); p++)
    {
      for (int i=0; i<2; i++)
	{
	  EdgeVector edge;
	  int hang   = p->second.hanging();
	  int rneigh = p->second.rneighbour();
	  pair<int,int> cp = QuadLaO.GetChildEdges(edge,p->first,hang,rneigh,i);

	  int cellindex = cp.first;
	  int edgeindex = cp.second;

	  Edge E(cellindex,edgeindex);

	  quads[cellindex].edge(edgeindex) = n;

	  edges[n] = E;

	  int bigslave = p->second.cneighbour();

	  if (bigslave>=0)
	    {
	      if (quads[bigslave].sleep())
		{
		  int slave = 0;
		  for (int j=0; j<4; j++)
		    {
		      slave = quads[bigslave].child(j);
		      edgeindex = QuadLaO.local_edge_index(slave,edge);
		      if (edgeindex>=0) break;
		    }
		  if (edgeindex<0)
		    {
		      cout << slave << " ###smallind 2 " << edgeindex << endl;
		      exit(1);
		    }
		  quads[slave].edge(edgeindex) = n;
		  edges[n].slave() = slave;
		  edges[n].LocalSlaveIndex() = edgeindex;
		}
	    }
	  n++;
	}
    }
}

/*---------------------------------------------------*/

void EdgeManager::OldHangings(HangContainer2d& hangset,
			      const IntSet& CellRefList)
{
  HangList<2>::iterator p = hangset.NotAnyMore().begin();

  for (;p!=hangset.NotAnyMore().end(); p++)
    {
      if (p->second.hanging()<0) continue;

      int bigslave  = p->second.cneighbour();
      int bigmaster = p->second.rneighbour();
      int hang      = p->second.hanging();

      if (hang<0)      continue;
      if (bigslave<0)  continue;
      if (bigmaster<0) continue;

      if (!quads[bigslave].sleep())
	{
	  cout << "BIGSLAVE sleeps" << endl;
	  exit(1);	      
	}
      for (int i=0; i<2; i++)
	{
	  EdgeVector edge;

	  int rneigh = p->second.rneighbour();
	  pair<int,int> cp = QuadLaO.GetChildEdges(edge,p->first,hang,rneigh,i);

	  int gedge  = quads[cp.first].edge(cp.second);
	  if (gedge<0)
	    {
	      swap(bigmaster,bigslave);
	      swap(p->second.rneighbour(),p->second.cneighbour());

	      rneigh = p->second.rneighbour();
	      cp = QuadLaO.GetChildEdges(edge,p->first,hang,rneigh,i);

	      gedge  = quads[cp.first].edge(cp.second);
	      if (gedge<0)
		{
		  cout << "###gedge negativ:" << bigmaster << " " << bigslave << endl;
		}
	    }
	  int ledge  = -1;
	  int slave  = -1;
	  int master = cp.first;

	  for (int j=0; j<4; j++)
	    {
	      slave = quads[bigslave].child(j);
	      ledge = QuadLaO.local_edge_index(slave,edge);
	      if (ledge>=0) break;
	    }
	  if (ledge<0)
	    {
	      cout << slave << " ###ledge " << ledge << endl;
	      exit(1);
	    }
	  Edge& E = edges[gedge];

	  if (E.master()!=master) 
	    {
	      if (E.slave()==master)
		{
		  E.master()=master;
		  E.LocalMasterIndex() = E.LocalSlaveIndex();
		}
	      else
		{
		  cout << endl << "bad master" << master << ": ";
		  cout << E.master() << "-" << E.slave()<< endl;
		  exit(1);
		}
	    }
	  quads[slave].edge(ledge) = gedge;
	  E.slave()           = slave;
	  E.LocalSlaveIndex() = ledge;
	}
    }
}

/*---------------------------------------------------*/

void EdgeManager::SwappedEdges()
{
  int n = 0;
  int m = 0;
  for (int i=0; i<quads.size(); i++)
    {
      const Quad& q = quads[i];
      if (q.sleep()) continue;
      for (int e=0; e<4; e++)
	{
	  if (q.edge(e)<0) m++;
	}
    }
  if (m!=SwappedEdge.size())
    {
      cout << "ERROR: Inconsistency in SwappedEdges2d" << endl;
      cout << m << " " << SwappedEdge.size() << endl;
      exit(1);
    }
  for (int i=0; i<quads.size(); i++)
    {
      Quad& q = quads[i];
      if (q.sleep()) continue;
      for (int e=0; e<4; e++)
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

/*---------------------------------------------------*/

fixarray<2,int> EdgeManager::ChildrenOfEdge(int e) const
{
  int s  = edges[e].slave();
  int is = edges[e].LocalSlaveIndex();

  if(s<0)   
    {
      cerr << "EdgeManager::ChildrenOfEdge()";
      cerr << "no slave\n";
      abort();
    }
  fixarray<2,int> f;
  for(int ii=0;ii<2;ii++)
    {
      int ic = quad(s).child(QuadLaO.ChildsOfEdge(is,ii));
      f[ii] = quad(ic).edge(QuadLaO.ChildEdge(is));
    }
  return f;
}

/*---------------------------------------------------*/

void EdgeManager::DeleteEdges()
{
  compress(edges,eo2n);

  for (int i=0; i<quads.size(); i++)
    {
      Quad& Q = quads[i];
      if (co2n[i]<0)
	{
	  Q.edges() = -1;
	}
      else
	{
	  for (int e=0; e<4; e++)
	    {
	      int ne = eo2n[ Q.edge(e) ];
	      if (ne<0)
		{
		  cout << "\neo2n " << ne;
		  exit(1);
		}
	      Q.edge(e) = ne;
	    }
	}
    }
}

/*---------------------------------------------------*/

bool EdgeManager::EdgeIsHanging(int e) const
{
  int m = edges[e].master();
  int s = edges[e].slave();
  if (s<0) return 0;
  if (quad(m).sleep() && !quad(s).sleep() ) return 1;
  if (quad(s).sleep() && !quad(m).sleep() ) return 1;
  return 0;
}

/*---------------------------------------------------*/

bool EdgeManager::EdgeIsHanging(const Edge& e) const
{
  int m = e.master();
  int s = e.slave();
  if (s<0) return 0;
  if (quad(m).sleep() && !quad(s).sleep() ) return 1;
  if (quad(s).sleep() && !quad(m).sleep() ) return 1;
  return 0;
}

/*---------------------------------------------------*/

void EdgeManager::SortHangings()
{
  // edges with hanging nodes swapped to the end of list

  vector<int> perm(edges.size());
  iota(perm.begin(),perm.end(),0);
  stable_sort(perm.begin(),perm.end(),HangEdgeSort(*this));
  stable_sort(edges.begin(),edges.end(),HangEdgeSort2(*this));
 
  int i=edges.size();
  while(EdgeIsHanging(i-1)) i--;

  vector<int> permi(perm.size());
  for(int i=0;i<perm.size();i++) permi[perm[i]] = i;

  for(int i=0;i<quads.size();i++)
    {
      for(int ii=0;ii<4;ii++)  quads[i].edge(ii) = permi[quad(i).edge(ii)];
    }
  for (int i=0; i<eo2n.size(); i++)
    {
      if (eo2n[i]>=0) eo2n[i] = permi[eo2n[i]]; 
    }

  // master of each edge is allways the coarser quad

  for(int i=0;i<edges.size();i++)
    {
      int m = edges[i].master();
      int s = edges[i].slave();
      if(s<0) continue;
      if (quad(m).sleep() && !quad(s).sleep())
	{
	  swap(edges[i].master(),edges[i].slave());
	  swap(edges[i].LocalMasterIndex(),edges[i].LocalSlaveIndex());
	}
    }
}

/*---------------------------------------------------*/

void EdgeManager::InitEdges()
{
  HANGMAP H;

  EdgeVector e;
  for (int i=0; i<quads.size(); i++)
    {
      for (int j=0; j<4; j++)
	{
	  e[0] = quads[i].vertex(j);
	  e[1] = quads[i].vertex((j+1)%4);

	  if (e[1]<e[0])  swap(e[0],e[1]);

	  HANGMAP::iterator yes = H.find(e);

	  if (yes!=H.end())
	    {
	      int k = yes->second;
	      edges[k].slave() = i;
	      edges[k].LocalSlaveIndex()=j;
	      quads[i].edge(j) = k;
	      H.erase(yes);
	    }
	  else
	    {
	      Edge E(i,j);
	      int n = edges.size();
	      edges.push_back(E);
	      H.insert(make_pair(e,n));
	      quads[i].edge(j) = n;
	    }
	}
    }
}

/*---------------------------------------------------*/

void EdgeManager::NeighbourTester() const
{
  int n = quads.size();
  vector<fixarray<4,int> >  vecino(n);

  for (int q=0; q<quads.size(); q++)
    {
      for (int e=0; e<4; e++)
	{
	  int ind = quads[q].edge(e);
	  if (ind<0)
	    {
	      cout << "* Q=" << q << " edge=" << e << " " << ind << endl;
	    }
	}
    }
  for (int q=0; q<quads.size(); q++)
    {
      for (int e=0; e<4; e++)
	{
	  int ind = quads[q].edge(e);
	  int m = edges[ind].master();
	  int s = edges[ind].slave();
	  if (m==q)
	    {
	      vecino[q][e]=s;
	    }
	  else if (s==q)
	    {
	      vecino[q][e]=m;
	    }
	  else
	    {
	      vecino[q][e]=-2;
	      cout << "*q=" << q << " el=" << e << " eg=" << ind << ": m=" << m << " s=" << s << endl;
	    }
	}
    }
}
}

/*------------------------------------------------------*/

#undef HANGMAP
