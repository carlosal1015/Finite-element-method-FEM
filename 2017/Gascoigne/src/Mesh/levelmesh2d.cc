/**
*
* Copyright (C) 2004, 2006, 2007, 2008 by the Gascoigne 3D authors
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


#include  "levelmesh2d.h"
#include  "levelsorter.h"
#include  "leveljumper.h"
#include  "set2vec.h"
#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_set>
#define HASHSET std::tr1::unordered_set
#else
#include  <ext/hash_set>
#define HASHSET __gnu_cxx::hash_set
#endif

using namespace std;
 
/*---------------------------------------------------*/

namespace Gascoigne
{
LevelMesh2d::LevelMesh2d(const HierarchicalMesh* hmp) 
  : Index()
{
  HMP = dynamic_cast<const HierarchicalMesh2d*>(hmp);
}

/*---------------------------------------------------*/

LevelMesh2d::~LevelMesh2d() 
{}

/*---------------------------------------------------*/

bool LevelMesh2d::EdgeIsHangingGlobalIndex(int i) const
{
  int igm = HMP->edge(i).master();
  int igs = HMP->edge(i).slave();

  // rand oder kleien hedges
  if(igs<0) return 0;

  bool m_in_lmesh = (Quadg2l().find(igm)!=Quadg2l().end());
  bool s_in_lmesh = (Quadg2l().find(igs)!=Quadg2l().end());

  if(m_in_lmesh && s_in_lmesh) return 0;

  int ivg = HMP->NodeOnEdge(i);
  if(Vertexg2l().find(ivg)==Vertexg2l().end()) return 0;

  return 1;
}

/*---------------------------------------------------*/

void LevelMesh2d::BasicInit(const IntSet& newq, const IntSet& oldq) 
{
  // doch sortiert

  int n = newq.size()+oldq.size();
  Index::Quadl2g().memory(n);
  IntVector::const_iterator p = set_union(newq.begin(),newq.end(),oldq.begin(),oldq.end(),
					  Index::Quadl2g().begin());
  n = p-Index::Quadl2g().begin();

  InitCells(n);
  InitNodes(n);
  //  InitEdges(n);
}

/*-----------------------------------------*/

void LevelMesh2d::InitCells(int n)
{
  Index::Quadl2g().memory(n);

  sort(Index::Quadl2g().begin(), Index::Quadl2g().end(), LevelSorter2d(*HMP));

  Index::InitQuads();
}

/*-----------------------------------------*/

void LevelMesh2d::InitNodes(int n)
{
  IntSet nodes;
  for(int i=0;i<n;i++)
    {
      int ind = Index::Quadl2g()[i];
      for(int ii=0;ii<4;ii++)  
	{
	  nodes.insert(HMP->vertex_of_cell(ind,ii));
	}
    }
  Index::InitNodes(nodes);
}

/*-----------------------------------------*/

void LevelMesh2d::InitEdges(int n)
{
  // edges
  IntSet edges;
  for(int i=0;i<n;i++)
    {
      const Quad& Q = HMP->quad(Index::Quadl2g()[i]);
      for(int ii=0;ii<4;ii++)  
	{
	  edges.insert(Q.edge(ii));
	}
    }
  Index::InitEdges(edges);

  // sortiere haengende edges nach hinten

  stable_sort(Edgel2g().begin(),Edgel2g().end(),HangEdgeSort3(*this));

  Edgeg2l().clear();
  for(int i=0;i<Edgel2g().size();i++)
    {
      Edgeg2l().insert(make_pair(Edgel2g()[i],i));
    }
}

/*-----------------------------------------*/

bool LevelMesh2d::BuildFathers(set<int>&  Vaeter) const
{
  for(int i=0; i<ncells();i++)
    {
      const Quad& q = quad(i);
      int findex = q.father();
      if(findex==-1) 
	{
	  return 0;
	}

      const Quad& qf = HMP->quad(findex);
      for(int ii=0;ii<qf.nchilds();ii++) {
	int cindex = qf.child(ii);
	if(Quadg2lCheck(cindex)==-2) {
	  return 0;
	}
      }
      Vaeter.insert(findex);
    }
  return 1;
}

/*-----------------------------------------*/

bool LevelMesh2d::ConstructCellIndOfPatch(IntVector& dst) const
{
  set<int>  Vaeter;
  BuildFathers(Vaeter);

  int nq = ncells()/4;
  dst.reservesize(nq);

  int count=0;
  set<int>::const_iterator pf = Vaeter.begin();
  while(pf!=Vaeter.end())
    {
      int findex = *pf;
      dst[count] = findex;
     
      count++;
      pf++;
    }
  return 1;
}


/*-----------------------------------------*/

void LevelMesh2d::ConstructIndOfPatch(nvector<IntVector>& dst) const
{
  set<int>  Vaeter;
  BuildFathers(Vaeter);

  int nq = ncells()/4;
  dst.reserve(nq);
  dst.resize (nq,IntVector(9));

  int count=0;
  set<int>::const_iterator pf = Vaeter.begin();
  while(pf!=Vaeter.end())
    {
      int findex = *pf;
      const Quad& qf = HMP->quad(findex);

      fixarray<4,int> FineQuads;
      for(int ii=0;ii<qf.nchilds();ii++) 
	{
	  FineQuads[ii] = qf.child(ii);
	}
      dst[count][0] = Vertexg2l( HMP->quad(FineQuads[0]).vertex(0) );
      dst[count][1] = Vertexg2l( HMP->quad(FineQuads[0]).vertex(1) );
      dst[count][2] = Vertexg2l( HMP->quad(FineQuads[1]).vertex(1) );
      dst[count][3] = Vertexg2l( HMP->quad(FineQuads[0]).vertex(3) );
      dst[count][4] = Vertexg2l( HMP->quad(FineQuads[0]).vertex(2) );
      dst[count][5] = Vertexg2l( HMP->quad(FineQuads[1]).vertex(2) );
      dst[count][6] = Vertexg2l( HMP->quad(FineQuads[3]).vertex(3) );
      dst[count][7] = Vertexg2l( HMP->quad(FineQuads[3]).vertex(2) );
      dst[count][8] = Vertexg2l( HMP->quad(FineQuads[2]).vertex(2) );
     
      count++;
      pf++;
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::ConstructHangingStructureQuadratic(QuadraticHNStructure3& hnq2) const
{
  hnq2.clear();
  set<int> habschon; 
  for (int i=0;i<ncells();i++)
    {
      const Quad& q = HMP->quad(Quadl2g(i));
      int father = q.father();
      if (father<0) continue;
      if (habschon.find(father)!=habschon.end()) continue;
      habschon.insert(father);
      for(int in=0; in<4; in++)
	{
	  int neighbour = HMP->neighbour(father,in);
	  if (neighbour<0) continue;
	  const Quad& qfn = HMP->quad(neighbour);
	  if (qfn.nchilds()==0) continue;

	  fixarray<2,int> childs;
	  int ne = in;
	  
	  {
	    int start = 0;
	    int neighbourneighbour = HMP->neighbour(neighbour,ne);
	    while ((neighbourneighbour!=father) && (start<4))
	      {
		start++;
		ne = (ne+1)%4;
		neighbourneighbour = HMP->neighbour(neighbour,ne);
	      }
	    
	    assert(neighbourneighbour==father);
	  }	  
	  HMP->QuadLawOrder().childs_of_edge(childs,qfn,ne);
	  int child = childs[0];
	  if (Quadg2lCheck(child)>=0) continue;
	  const Quad& qfc = HMP->quad(child);
	  
	  if (qfc.nchilds()==0) continue;

	  int enkel = qfc.child(0);
	  if (Quadg2lCheck(enkel)<0) continue;
	  
	  // jetzt haengt er
	  int hn = Vertexg2l( HMP->QuadLawOrder().edge_vertex(qfc,ne) );
	  fixarray<3,int> F;
	  F[0] = qfn[ne];
	  F[1] = HMP->QuadLawOrder().edge_vertex(qfn,ne);
	  F[2] = qfn[(ne+1)%4];
	  
	  assert((qfc[0]==F[0]) || (qfc[1]==F[0]) ||
		 (qfc[2]==F[0]) || (qfc[3]==F[0]));
	  
	  for (int k=0; k<3; k++)  F[k] = Vertexg2l(F[k]);
	  
	  hnq2.insert(make_pair(hn,F));
	  
	  const Quad& qfc2 = HMP->quad(childs[1]);
	  hn = Vertexg2l(HMP->QuadLawOrder().edge_vertex(qfc2,ne));
	  
	  swap(F[0],F[2]);
	  
	  hnq2.insert(make_pair(hn,F));
	}
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::ConstructHangingStructureQuartic(QuarticHNStructure5& hnq4) const
{
  hnq4.clear();
  assert(HMP->patchdepth()>=2);
  int count = 0;
  std::set<int> habschon;
  for (int i=0;i<ncells();++i)
    {
      const Quad& q_c = HMP->quad(Quadl2g(i));
      int father_c    = q_c.father();
      assert(father_c>=0);
      const Quad& f_c = HMP->quad(father_c);
      int opa_c       = f_c.father();
      assert(opa_c>=0);
      

		// wurde der opa schon behandelt?
      if (habschon.find(opa_c)!=habschon.end()) continue;
      habschon.insert(opa_c);
      
		// alle Nachbarn
      for (int in=0;in<4;++in)
	{
	  int opa_r   = HMP->neighbour(opa_c,in);
		    // Nachbar existiert nicht
	  if (opa_r<0) continue;
	  const Quad& o_r = HMP->quad(opa_r);
		    // muss einmal verfeinert sein
	  assert(o_r.nchilds()!=0);

		    // nachbar-nachbar finden
	  int ne = in;
	  {
	    int start =0;
	    int neighbourneighbour = HMP->neighbour(opa_r,ne);
	    while ((neighbourneighbour!=opa_c) && (start<4))
	      {
		start++;
		ne = (ne+1)%4;
		neighbourneighbour = HMP->neighbour(opa_r,ne);
	      }
	    assert(neighbourneighbour==opa_c);
	  }
		    // beide Vaeter entlang der gemeinsamen Kante
	  fixarray<2,int> fathers_r;
	  HMP->QuadLawOrder().childs_of_edge(fathers_r,o_r,ne);

		    // wenn gleiches Level, dann haengt nix.
	  const Quad& f0_r = HMP->quad(fathers_r[0]);
	  const Quad& f1_r = HMP->quad(fathers_r[1]);
		    // groeber?
	  if (f0_r.nchilds()==0) continue;
		    // gleiches level?
	  const Quad& q0_r = HMP->quad(f0_r.child(0));
	  if (q0_r.nchilds()==0) continue;

		    // der Enkel muss aktiv sein.
	  assert(Quadg2lCheck(q0_r.child(0))>=0);

		    // Jetzt haengt ne Menge

		    // die vier enkel verfeinert
	  vector<fixarray<2,int> > enkels_r(2);
	  HMP->QuadLawOrder().childs_of_edge(enkels_r[0],f0_r,ne);
	  HMP->QuadLawOrder().childs_of_edge(enkels_r[1],f1_r,ne);
	  assert(enkels_r[0][0]>=0);assert(enkels_r[0][1]>=0);assert(enkels_r[1][0]>=0);assert(enkels_r[1][1]>=0);
	  
 
		    // die 5 Knoten
	  fixarray<6,int> tmp;
	  tmp[0]=HMP->quad(enkels_r[0][0])[ne];
	  tmp[1]=HMP->quad(enkels_r[0][1])[ne];
	  tmp[2]=HMP->quad(enkels_r[1][0])[ne];
	  tmp[3]=HMP->quad(enkels_r[1][1])[ne];
	  tmp[4]=HMP->quad(enkels_r[1][1])[(ne+1)%4];
	  

		    // die haengen
	  fixarray<4,int> hn;
	  hn[0] = HMP->QuadLawOrder().edge_vertex(HMP->quad(enkels_r[0][0]),ne);
	  hn[1] = HMP->QuadLawOrder().edge_vertex(HMP->quad(enkels_r[0][1]),ne);
	  hn[2] = HMP->QuadLawOrder().edge_vertex(HMP->quad(enkels_r[1][0]),ne);
	  hn[3] = HMP->QuadLawOrder().edge_vertex(HMP->quad(enkels_r[1][1]),ne);

	  tmp[5]=0;   hnq4.insert(std::make_pair(hn[0],tmp));
	  tmp[5]=1;   hnq4.insert(std::make_pair(hn[1],tmp));
	  std::swap(tmp[0],tmp[4]);
	  std::swap(tmp[1],tmp[3]);
	  tmp[5]=1;   hnq4.insert(std::make_pair(hn[2],tmp));
	  tmp[5]=0;   hnq4.insert(std::make_pair(hn[3],tmp));
	  

	  ++count;
	  
	}
      
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::check_leveljump() const
{
  LevelJumper  Phi;
  for(int c=0;c<ncells();c++)
    {
      Phi.update(quad(c));
    }
  assert(! Phi.check());
  //if(Phi.check()) cerr << "LevelMesh2d::check_leveljump() aenderb\n";
}

/*---------------------------------------------------*/

void LevelMesh2d::fill_opis(IntSet& dst, IntSet& oldquads) const
{
  for(int i=0; i<ncells(); i++)
    {
      const Quad& Q = quad(i);

      int f = Q.father();
      assert(f>=0);

      int opi = HMP->quad(f).father();

      if (opi<0) 
	{
	  int j = Quadl2g(i);
	  oldquads.insert(j);
	}
      else
	{
	  dst.insert(opi);
	}
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::fill_childs(IntSet& dst, const Quad& Q) const
{
  for (int ii=0; ii<Q.nchilds(); ii++)
    {
      int qccindex = Q.child(ii);
      dst.insert(qccindex);
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::fill_enkel(IntSet& oldquads, const Quad& Q) const
{
  for (int ii=0; ii<Q.nchilds(); ii++)
    {
      int qccindex = Q.child(ii);
      const Quad& qcc = HMP->quad(qccindex);
      for(int iii=0;iii<qcc.nchilds();iii++)
	{
	  int qcindex = qcc.child(iii);
	  if(Quadg2lCheck(qcindex)>=0)
	    {
	      oldquads.insert(qcindex);
	    }
	}
    }
}

/*---------------------------------------------------*/

bool LevelMesh2d::EnkelUniform(const Quad& Q) const
{
  for (int ii=0; ii<Q.nchilds(); ii++)
    {
      int qccindex = Q.child(ii);
      const Quad& qcc = HMP->quad(qccindex);
      for(int iii=0;iii<qcc.nchilds();iii++)
	{
	  int qcindex = qcc.child(iii);
	  if(Quadg2lCheck(qcindex)==-2)
	    {
	      return 0;
	    }
	}
    }
  return 1;
}

/*---------------------------------------------------*/

void LevelMesh2d::construct_lists(IntSet& newquads, IntSet& oldquads) const
{
  newquads.clear();  
  oldquads.clear();

  check_leveljump();

  set<int>     Opis;
  fill_opis(Opis,oldquads);
  for(set<int>::const_iterator p=Opis.begin();p!=Opis.end();p++)
    {
      const Quad& Q = HMP->quad(*p);
      
      if ( EnkelUniform(Q) )
	{
	  fill_childs(newquads,Q);
	}
      else
	{
	  fill_enkel(oldquads,Q);
	}
    }

  // Iteration zum Regulaer machen (LevelJump)
  while(1)
    {
      LevelJumper  Phi;
      set<int>::const_iterator p;
      for(p=newquads.begin(); p!=newquads.end(); p++)
	{
	  Phi.update(HMP->quad(*p));
	}
      for(p=oldquads.begin(); p!=oldquads.end(); p++)
	{
	  Phi.update(HMP->quad(*p));
	}
      if(!Phi.check()) break;

      //int rep=0;
      IntSet help(newquads);
      for(p=newquads.begin(); p!=newquads.end(); p++)
	{
	  const Quad& q = HMP->quad(*p);
	  if (!Phi.VertexOK(q))
	    {
	      //rep++;
	      const Quad& qf = HMP->quad(q.father());
	      for(int ii=0;ii<4;ii++)
		{
		  help.erase(qf.child(ii));
		}
	      fill_enkel(oldquads,qf);
	    }
	}
      newquads = help;
      //cerr << "\t Regular Iteration\t" << count++ << " " << rep << endl;
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::InitBoundaryHandler(BoundaryIndexHandler& BI,const PatchIndexHandler& PIH) const
{
  IntSet blines;
  for(int i=0;i<HMP->nblines();i++)
    {
      int q = HMP->bline(i).of_quad();
      if(Quadg2l().find(q)!=Quadg2l().end())
	{
	  blines.insert(i);
	}
    }

  // which colors are there ?
  BI.clear();
  for(IntSet::const_iterator p=blines.begin();
      p!=blines.end(); p++)
    {
      const BoundaryLine& bl = HMP->bline(*p);
      int col = bl.material();

      BI.GetColors().insert(col);
    }
  IntVector colorvec;
  Set2Vec(colorvec,BI.GetColors());

  // compute inverse positions
  map<int,int> inv;

  for (int i=0; i<colorvec.size(); i++)
    {
      inv.insert(make_pair(colorvec[i],i));
    }
  
  int nc = colorvec.size(); 
  vector<set<int>  > H1(nc);  // for verteces
  // for cells and local indices
  vector<set<fixarray<2,int> > >  H2(nc); 

  for(IntSet::const_iterator q=blines.begin();
      q!=blines.end(); q++)
    {
      const BoundaryLine& bl = HMP->bline(*q);
      int col = bl.material();

      map<int,int>::const_iterator p = inv.find(col);
      if (p==inv.end())
	{
	  cout << "LevelMesh2d::BuildVertexOnBoundary()"<< endl;
	  abort();
	}
      int pos = p->second;

      for(int ii=0;ii<2;ii++)
	{
	  int vindex = Vertexg2l(bl.vertex(ii));
	  H1[pos].insert(vindex);
	}
      fixarray<2,int> ind;
      ind[0] = Quadg2l(bl.of_quad());
      ind[1] = bl.edge_in_quad();
      H2[pos].insert(ind);
    }
  BI.CopySetToVector(H1,colorvec,BI.GetVertex());

  for (int i=0; i<H2.size(); i++)
    {
      IntVector v1(H2[i].size());
      IntVector v2(H2[i].size());
      int j = 0;
      
      set<fixarray<2,int> >::const_iterator p;
      for (p=H2[i].begin(); p!=H2[i].end(); p++)
	{
	  v1[j] = (*p)[0];
	  v2[j] = (*p)[1];
	  j++;
	}
      int color = colorvec[i];

      BI.GetCell().insert(make_pair(color,v1));
      BI.GetLocal().insert(make_pair(color,v2));
    }
  
  const nvector<IntVector>& patch2cell = PIH.GetAllPatch2Cell();
  
  nvector<int> cell2patch(PIH.npatches()<<2);
  for (int p=0;p<patch2cell.size();++p)
    for (int i=0;i<patch2cell[p].size();++i)
      cell2patch[patch2cell[p][i]]=p;
 
  for(IntSet::const_iterator c=BI.GetColors().begin(); c != BI.GetColors().end(); c++)
  {
    int col = *c;
    const IntVector& cells = BI.Cells(col);
    const IntVector& locals= BI.Localind(col);
    HASHSET<int> habschon;

    IntVector p1;
    IntVector p2;
      
    for(int i=0; i < cells.size(); i++)
    {
      int iq = cells[i];
      int ip = cell2patch[iq];
      int ile = locals[i];

      //gabs den patch schon
      if(habschon.find((ip<<2)+ile) != habschon.end()) continue;
      habschon.insert((ip<<2)+ile);
      p1.push_back(ip);
      p2.push_back(ile);
    }
    BI.GetPatch().insert(make_pair(col,p1));
    BI.GetLocalPatch().insert(make_pair(col,p2));
  }

}
}
#undef HASHSET
/*--------------------------------------------------------------*/

