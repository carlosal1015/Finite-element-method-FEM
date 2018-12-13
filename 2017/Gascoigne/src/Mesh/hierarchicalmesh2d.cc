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


#include  "hierarchicalmesh2d.h"
#include  "vecalgo.h"
#include  "deletecells.h"
#include  "set2vec.h"
#include  "levelcomparer2d.h"
#include  "edgemanager.h"
#include  "stlio.h"
#include  "regular_update.h"
#include  "coarsehierarchicalmesh2d.h"

#include  <set>
#include  <iterator>
#include  <fstream>
#include  <stdio.h>
#include  "giota.h"

using namespace std;

/*------------------------------------------------------*/

namespace Gascoigne
{
HierarchicalMesh2d::HierarchicalMesh2d() 
  : HierarchicalMesh(), QuadLaO(quads) {}

/*------------------------------------------------------*/

HierarchicalMesh2d::HierarchicalMesh2d(const ParamFile* paramfile) 
  : HierarchicalMesh(), QuadLaO(quads)
{
  BasicInit(paramfile);
}

/*------------------------------------------------------*/

HierarchicalMesh2d::HierarchicalMesh2d(const HierarchicalMesh2d& H)
  : HierarchicalMesh(), QuadLaO(quads)
{
  *this = H;
}

/*------------------------------------------------------*/

HierarchicalMesh2d& HierarchicalMesh2d::operator=(const HierarchicalMesh2d& H)
{
  HierarchicalMesh::operator= (H);
  // copy all data
  vertexs2d = H.vertex2d();
  quads     = H.quad();
  Blines    = H.bline();
  LineHang  = H.linehang();
  // line shapes werden nicht kopiert sonst geht patchrefine in die britze
  quadofcurved = H.GetQuadOfCurved();

  return *this;
}

/*------------------------------------------------------*/

pair<int,int> HierarchicalMesh2d::GetBoundaryInformation(int i) const
{
  int material = -1; 
  int le = -1;
  int ib = GetBoundaryCellOfCurved(i);
  if (ib>=0)
    {
      material = bline(ib).material();
      le       = bline(ib).edge_in_quad();
    }
  return make_pair(material,le);
}

/*------------------------------------------------------*/

int HierarchicalMesh2d::FindPatchDepth() const
{
  //simple version, sucht nur p=1, p=0
  for(int i=0;i<ncells();i++)
    {
      const Quad& q = quad(i);
      if(q.sleep()) continue;
      int father = q.father();
      if(father==-1) return 0;
      const Quad& qf = quad(father);
      for(int ii=0;ii<4;ii++) 
	{
	  if(quad(qf.child(ii)).sleep()) return 0;
	}
    }
  return 1;
}

/*------------------------------------------------------*/

const BoundaryFunction<2>* HierarchicalMesh2d::line_shape(int i) const
{
  if (GetCurvedShapes().empty()) return NULL;
  if (GetCurvedShapes().Curved(i)) return &(GetCurvedShapes().GetShape(i));
  return NULL;
}

/*------------------------------------------------------*/

set<int> HierarchicalMesh2d::GetColors() const
{
  set<int> coleur;
  for(int i=0;i<nblines();i++)
    {
      coleur.insert(bline(i).material());
    }
  return coleur;
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::InitQuadOfCurved()
{
  quadofcurved.clear();

  if (GetCurvedShapes().empty()) return;

  for(int il=0;il<nblines();++il)
    {
      const BoundaryLine& B = bline(il);

      if (GetCurvedShapes().Curved(B.material()))
	{
	  int iq = B.of_quad();
	  quadofcurved.insert(make_pair(iq,il));
	}
    }
}


/*------------------------------------------------------*/

int    HierarchicalMesh2d::Vater(const int i) const
{
  return quad(i).father();
}


/*------------------------------------------------------*/

IntVector HierarchicalMesh2d::Nachkommen(const int i) const
{
  IntVector k = Kinder(i);
  if (k.size()==0) return k;
  IntVector k1;
  int ks=k.size();
  for (int i=0;i<ks;++i)
    {
      IntVector k1 = Nachkommen(k[i]);
      for (int j=0;j<k1.size();++j) k.push_back(k1[j]);
    }
  return k;
}

/*------------------------------------------------------*/

int HierarchicalMesh2d::nactivedescendants(int i)      const
{
  if (!quads[i].sleep())
    return 1;
  int k=0;
  for (int j=0;j<quads[i].nchilds();++j)
    k+=nactivedescendants(quads[i].child(j));
  return k;
}

/*------------------------------------------------------*/

IntVector HierarchicalMesh2d::GetVertices(int c) const
{  
  IntVector v;
  for (int i=0;i<quads[c].nvertexs();++i)
    v.push_back(quads[c][i]);
  return v;
}

/*------------------------------------------------------*/

IntVector HierarchicalMesh2d::Kinder     (const int i) const
{
  IntVector k = quad(i).childs();
  return k;
} 

/*------------------------------------------------------*/


IntVector HierarchicalMesh2d::Geschwister(const int i) const
{
  const Quad& q = quad(i);
  int father = q.father();
  if (father==-1)
    {
      IntVector n(1,i);
      return n;
    }
  return Kinder(father);
}


/*------------------------------------------------------*/

fixarray<2,int> HierarchicalMesh2d::ChildrenOfEdge(int e) const
{
  int s  = edge(e).slave();
  int is = edge(e).LocalSlaveIndex();

  assert(s>=0) ;  

  fixarray<2,int> f;
  for(int ii=0;ii<2;ii++)
    {
      int ic = quad(s).child(QuadLaO.ChildsOfEdge(is,ii));
      f[ii] = quad(ic).edge(QuadLaO.ChildEdge(is));
    }

  return f;
}

/*---------------------------------------------------*/

int  HierarchicalMesh2d::NodeOnEdge(int e) const
{
  // only for real hanging nodes
  int m  = edge(e).master();
  int im = edge(e).LocalMasterIndex();
  int s  = edge(e).slave();
  int is = edge(e).LocalSlaveIndex();

  if(s<0) 
    {
      assert(quad(im).sleep());

      return QuadLaO.edge_vertex(quad(m),im);
    }

  return QuadLaO.edge_vertex(quad(s),is);
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::ghostglobalcoarse(HangContainer2d& hangset,
					 const IntSet& cellcoarse)
{
  EdgeVector lineglob;

  IntVector coarse;
  Set2Vec(coarse,cellcoarse);
  LevelComparer2d  lc(*this,coarse);
  
  IntVector BTS(lc.size()); iota(BTS.begin(),BTS.end(),0);
  sort(BTS.begin(),BTS.end(),CompareObjectBigToSmall<LevelComparer2d> (lc));
  
  for (int cp=0; cp<coarse.size(); cp++)
    {
      int f = coarse[BTS[cp]];
      for(int edge=0;edge<4;++edge)
	{
	  QuadLaO.global_edge_unsorted(lineglob,quad(f),edge);

	  int ve = QuadLaO.edge_vertex(quad(f),edge);
	  hangset.ghost_coarse(lineglob,f,ve);
	}
    }
  ghost_fill_neighbours2d();
  hangset.make_consistent();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::ghost2d(HangContainer2d& hangset,
			     const IntSet& cellref, 
			     const IntSet& cellcoarse)
{
  EdgeVector lineglob;
   
  for (IntSetIt cp=cellcoarse.begin(); cp!=cellcoarse.end(); ++cp)
    {
      int f = *cp;
      for(int edge=0;edge<4;++edge)
	{
	  QuadLaO.global_edge_unsorted(lineglob,quad(f),edge);

	  int ve = QuadLaO.edge_vertex(quad(f),edge);
	  hangset.ghost_coarse(lineglob,f,ve);
	}
    }
  for (IntSetIt cp=cellref.begin(); cp!=cellref.end(); ++cp)
    {
      int f = *cp;
      for(int edge=0;edge<4;++edge)
	{
	  QuadLaO.global_edge_unsorted(lineglob,quad(f),edge);

	  hangset.ghost_refine(lineglob,f);
	}
    }
  ghost_fill_neighbours2d();
  hangset.make_consistent();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::prepare2d(const IntVector& cell_ref, 
			       const IntVector& cell_coarse,
			       IntSet& CellRefList, IntSet& CellCoarseList)
{
  /* copies cell_ref into CellRefList without duplets  */

  IntVector::const_iterator  cp=cell_ref.begin();
  while(cp!=cell_ref.end())
    {
      int c = *cp;
      if ( (c>=0) && (c<quads.size()) )
	{
	  if(!quads[c].sleep())
	    {
	      CellRefList.insert(c);
	    }
	}
      ++cp;
    }

  /* copies cell_coarse into CellCoarseList without duplets  
     checks if :  -- coarse cell in refine list
                  -- coarse cell is cneighbour of a hang */

  IntSet  help;  

  for(IntVector::const_iterator cp = cell_coarse.begin(); 
      cp!=cell_coarse.end(); ++cp)
    {
      int ic = *cp;
      if((ic<0)||(ic>=quads.size())) continue;

      Quad& q = quads[ic];
      if (  q.sleep() )                              continue;
      if ( !q.level() )                              continue;
      if ( CellRefList.find(ic)!=CellRefList.end() ) continue;

      help.insert(ic);
    }

  LineHangList::const_iterator Lp;

  for(Lp = LineHang.begin(); Lp!=LineHang.end();Lp++)
    {
      int cn = Lp->second.cneighbour();
      if( help.find(cn)!=help.end() )    help.erase(cn);
    }

  multiset<int> ff;
  
  for (IntSetCIt hp = help.begin();
       hp != help.end(); ++hp)
    {
      ff.insert(quad(*hp).father());
    }

  for (multiset<int>::iterator fp = ff.begin();
       fp != ff.end(); ++fp)
    {
      if (ff.count(*fp)==4) 
	{
	  CellCoarseList.insert(*fp);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::ghost_fill_neighbours2d()
{
  EdgeVector lineglob;
  for(int ic=0;ic<quads.size();ic++)
    {
      for(int i=0;i<4;i++)
	{
	  QuadLaO.global_edge_unsorted(lineglob,quad(ic),i);
	  LineHangList::iterator p = LineHang.find(lineglob);

	  if(p!=LineHang.end())
	    {
	      int cn = p->second.cneighbour();
	      int rn = p->second.rneighbour();

	      // coarse neighbour

	      if( (cn==-1) && (rn!=ic) )
		{
		  p->second.cneighbour() = ic;
		}
	      
	      // refined neighbour

	      if( (rn==-1) && (cn!=ic) )
		{
		  p->second.rneighbour() = ic;
		}
	    }
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::basic_fill_neighbours2d()
{
  EdgeVector lineglob;
  for(int ic=0;ic<quads.size();ic++)
    {
      for(int i=0;i<4;i++)
	{
	  QuadLaO.global_edge_unsorted(lineglob,quad(ic),i);
	  LineHangList::iterator p = LineHang.find(lineglob);
	  
	  if(p!=LineHang.end())
	    {
	      int cn = p->second.cneighbour();
	      int rn = p->second.rneighbour();

	      if( (cn==-1) && (!quads[ic].sleep()) )
		{
		  p->second.cneighbour() = ic;
		}
	      else if(rn==-1)
		{
		  p->second.rneighbour() = ic;
		}
	    }
	 }
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::basic_refine2d(HangContainer2d& hangset,
				      const IntSet& CellRefList,
				      const IntSet& CellCoarseList)
{
  int ov = nnodes();
  int oc = quads.size();
  int oe = edges.size();

  int csub = 4*CellCoarseList.size();
  int cadd = 4*CellRefList.size();

  int vsub = hangset.NToBeDeleted()+CellCoarseList.size();
  int vadd = hangset.NToBeCreated()+CellRefList.size();

  int cdiff = cadd - csub;
  int vdiff = vadd - vsub;

  int nv = ov + vdiff;
  int nc = oc + cdiff;

//   cerr << "nv ov vdiff " << nv <<" "<< ov<<" "<< vdiff<<endl;

  clear_transfer_lists();

  IntVector cdel,vdel,edel;

  for(IntSetCIt p=CellCoarseList.begin();
      p!=CellCoarseList.end();p++)
    {
      const Quad& q = quad(*p);
      vdel.push_back(QuadLaO.middle_vertex(q));
      cdel.insert(cdel.end(),q.childs().begin(),q.childs().end());
    }
  hangset.load_elimination(vdel);
  hangset.NeighbourSwapper();

  //////
  EdgeManager EM(edges,quads,co2n,eo2n);
  //////

  EM.LoadEdgeElimination(edel,CellCoarseList,hangset);
 
  IntSet  LineRefList, LineCoarseList, ccdel;

  boundary_prepare2d(LineRefList, LineCoarseList, ccdel, hangset);
 
  transfer(oc,co2n,cdel);
  transfer(ov,vo2n,vdel);
  transfer(oe,eo2n,edel);

  IntVector cnew(cadd), vnew(vadd);
  iota(cnew.begin(),cnew.end(),oc-csub);
  iota(vnew.begin(),vnew.end(),ov-vsub);

  delete_vertexs2d  (vo2n);
  /////
  EM.DeleteEdges();
  /////
  delete_cells<Quad>(CellCoarseList,quads,co2n,vo2n);

  vertexs2d.reserve(nv);
  vertexs2d.resize(nv);
  quads    .reserve(nc);
  quads    .resize(nc);

  hangset.update_olds (vo2n,co2n);
  hangset.update_news (vnew, CellRefList.size());

  new_vertexs2d(hangset,vnew,CellRefList);
  new_quads    (hangset,cnew,vnew,ov,CellRefList);

  basic_fill_neighbours2d();

  new_boundary2d(LineRefList, LineCoarseList, ccdel);
  
  boundary_newton2d();
  InitQuadOfCurved();
  inner_vertex_newton2d(vnew,CellRefList);

  EM.Build(CellRefList,hangset);

}

/*---------------------------------------------------*/

void HierarchicalMesh2d::boundary_prepare2d(IntSet& LineRefList, 
					  IntSet& LineCoarseList,
					  IntSet& ccdel,
					  const HangContainer2d& hangset)
{
  EdgeVector   lineglob;
  for(int i=0; i<Blines.size(); i++)
    {
      const BoundaryLine& bl = Blines[i];
      
      lineglob[0] = bl.vertex(0);
      lineglob[1] = bl.vertex(1);
      sort(lineglob.begin(),lineglob.end());
      
      if (bl.sleep())
	{
	  if (hangset.ToBeDeleted(lineglob))
	    {
	      LineCoarseList.insert(i);
	      for (int j=0; j<2; j++) ccdel.insert(bl.child(j));
	    }
	}
      else
	{
	  if (hangset.ToBeCreated(lineglob))
	    {
	      LineRefList.insert(i);
	    }
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::new_boundary2d(IntSet& ref, IntSet& coarse,
				      IntSet& ccdel)
{
  int oc   = Blines.size();
  int csub = 2 * coarse.size();
  int cadd = 2 * ref.size();
  int nc   = oc+cadd-csub;

  IntVector lo2n;
  transfer(oc,lo2n,ccdel);
  delete_cells<BoundaryLine>(coarse,Blines,lo2n,vo2n);

  update_boundary_data2d(coarse);

  IntVector cnew(cadd);
  iota(cnew.begin(),cnew.end(),oc-csub);
  
  Blines.reserve(nc);
  Blines.resize(nc);

  new_lines(lo2n,cnew,ref);
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::new_quads(const HangContainer2d& hangset,
				 const IntVector& cnew, 
				 const IntVector& vnew, int nvold,
				 const IntSet& CellRefList)
{
  //cerr << "new_quads()" << endl;
  // neue zellen erzeugen
  // eintragen der "Vater-Vertexs" in den kindern

  int nci = 0;
  int ivm = 0;

  IntSetIt  cp=CellRefList.begin();
  EdgeVector lineglob;
  while(cp!=CellRefList.end())
    {
      int father = co2n[*cp++];

      assert(father>=0);

      vector<int>&  qc = quads[father].childs();
      qc.resize(4);

      int childlevel = quads[father].level()+1;
      for(int ic=0;ic<4;ic++)
	{
	  int inold = cnew[nci+ic];
	  qc[ic] = inold;
	  quads[inold].level()  = childlevel;
	  quads[inold].father() = father;
	  quads[inold].childs().resize(0);
	  quads[inold].edges() = -1;
	}
      QuadLaO.fill_corner_vertex_in_childs(quads[father]);
      QuadLaO.fill_middle_vertex_in_childs(quads[father],vnew[ivm]);
      ivm++;
      nci += 4;

      // Edge Vertex -- linehanginfo schon ok (hanging) !
      int ive(-1);
      for (int i=0; i<4; i++)
	{
	  QuadLaO.global_edge_unsorted(lineglob,quad(father),i);
	  
	  ive = hangset.vertex_index(lineglob);

	  assert(ive>=0);

	  QuadLaO.fill_edge_vertex_in_childs(quads[father],i,ive);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::inner_vertex_newton2d(const IntVector& vnew,
					     const IntSet& CellRefList)
{
  if (GetCurvedShapes().empty()) return;

  fixarray<4,int> v;
  fixarray<2,int> w; 
  IntSetIt  cp=CellRefList.begin();

  for (int i=0; i<CellRefList.size(); i++)
    {
      int oldind = *cp;
      int ind = co2n[*cp];
      cp++;
      if (quadofcurved.find(oldind)==quadofcurved.end()) continue;

      const Quad& f = quad(ind);

      // alter version
//       for (int j=0; j<4; j++) v[j] = QuadLawOrder().edge_vertex(f,j);
//       new_face_vertex2d(vnew[i],v);

      // neue Version
      //
      int bl = quadofcurved.find(oldind)->second;
      int ei = bline(bl).edge_in_quad();
      int ei2 = (ei+2)%4;
      w[0] = QuadLawOrder().edge_vertex(f,ei); 
      w[1] = QuadLawOrder().edge_vertex(f,ei2);
      new_edge_vertex2d(vnew[i],w);
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::new_vertexs2d(HangContainer2d& hangset, 
				     const IntVector& vnew,
				     const IntSet& CellRefList)
{
  int nv1 = CellRefList.size();

  IntSetIt  cp=CellRefList.begin();
  for (int i=0; i<nv1; i++)
    {
      int f  = co2n[*cp++];

      new_face_vertex2d(vnew[i],quads[f]);
    }
  HangList<2>::const_iterator p = hangset.Creating().begin();

  for (;p!=hangset.Creating().end(); p++)
    {
      new_edge_vertex2d(p->second.hanging(),p->first);
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::init_line(BoundaryLine& newline)
{
  fixarray<2,int> v,w;
  for (int i=0; i<quads.size(); i++)
    {
      for (int edge=0; edge<4; edge++)
	{
	  v[0] = quad(i).vertex(edge);
	  v[1] = quad(i).vertex((edge+1)%4);
	  if (newline==v)
	    {
	      newline.of_quad()      = i;
	      newline.edge_in_quad() = edge;
	      return;
	    }
	  w[0] = v[1];
	  w[1] = v[0];
	  if (newline==w)
	    {
	      newline.of_quad()      = i;
	      newline.edge_in_quad() = edge;
	      newline[0] = v[0];
	      newline[1] = v[1];
	      return;
	    }
	}
    }
  cerr << "Sophie im Brunnen !" << endl;
  cerr << newline;
  abort();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::update_boundary_data2d(const IntSet& LCoarse)
{
  int no = Blines.size() - 2*LCoarse.size();
  for(int i=0; i<no; ++i)
    {
      int oq = Blines[i].of_quad();
      assert(co2n[oq]>=0);

      Blines[i].of_quad() = co2n[oq];
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::boundary_newton2d()
{
  if (GetCurvedShapes().empty()) return;

  for (int i=0; i<Blines.size(); i++)
    {
      BoundaryLine& bl = Blines[i];
      int color = bl.material();
      
      if (GetCurvedShapes().Curved(color))
	{
	  GetCurvedShapes().newton(color, vertexs2d[bl.vertex(0)]);
	  GetCurvedShapes().newton(color, vertexs2d[bl.vertex(1)]);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::new_lines(const IntVector& lo2n, 
				 const IntVector& cnew,
				 const IntSet& LRef) 
{
  int nci = 0;

  for(IntSetCIt cp=LRef.begin(); cp!=LRef.end(); ++cp)
    {
      int father = lo2n[*cp];

      assert(father>=0);

      BoundaryLine& blf = Blines[father];
      // change father boundary
      vector<int>& qc = blf.childs();
      qc.resize(2);

      int edge = blf.edge_in_quad();
      int iq   = blf.of_quad();
      int vm   = QuadLaO.edge_vertex(quad(iq),edge);
      fixarray<2,int> chvec;
      QuadLaO.childs_of_edge(chvec,quad(iq),edge);

      for(int ic=0;ic<2;ic++)
	{
	  int inold = cnew[nci+ic];
	  // set childs in father
	  qc[ic]    = inold;
	  // set properties of childs
	  BoundaryLine& bl  = Blines[inold];
	  bl.level   () = blf.level()+1;
	  bl.father  () = father;
	  bl.material() = blf.material();
	  bl.childs().resize(0);

	  bl.vertex(ic)       =  blf.vertex(ic);
	  bl.vertex((ic+1)%2) =  vm;

	  bl.of_quad()      = chvec[ic];
	  bl.edge_in_quad() = QuadLaO.ChildEdge(edge);
	}

      nci += 2;
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::refine
(const IntVector& cell_ref, const IntVector& cell_coarse)
{
  IntSet CellRefList, CellCoarseList;
  _refine2d(CellRefList,CellCoarseList,cell_ref,cell_coarse);
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::_refine2d
(IntSet& CellRefList, IntSet& CellCoarseList, 
 const IntVector& cell_ref, const IntVector& cell_coarse)
{
  prepare2d(cell_ref,cell_coarse, CellRefList, CellCoarseList);

  while (regular_grid2d_three_refine(CellRefList)) {}
  while (regular_grid2d_three_coarse(CellRefList,CellCoarseList)) {}

  HangContainer2d hangset(LineHang);

  ghost2d(hangset,CellRefList,CellCoarseList);

  basic_refine2d(hangset,CellRefList,CellCoarseList);

  post_refine2d();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::delete_vertexs2d(const IntVector& vo2n)
{
  for(unsigned oi=0; oi<vo2n.size(); ++oi)
    {
      int ni = vo2n[oi];
      if(ni>=0)  
	{
	  vertexs2d[ni] = vertexs2d[oi];
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::post_refine2d()
{
  check_mesh2d();
  mnlevels = 0;
  for(int i=0;i<quads.size();i++)
    {
      mnlevels = Gascoigne::max_int(mnlevels,quads[i].level());
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::new_edge_vertex2d(int nv, const EdgeVector& v)
{
  vertexs2d[nv].equ(0.5, vertexs2d[v[0]], 0.5, vertexs2d[v[1]]);
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::new_face_vertex2d(int nv, const FaceVector& v)
{
  vertexs2d[nv].equ(0.25, vertexs2d[v[0]], 0.25, vertexs2d[v[1]],
		    0.25, vertexs2d[v[2]], 0.25, vertexs2d[v[3]]);
}

/*---------------------------------------------------*/

ostream& operator<<(ostream& os, const pair<fixarray<2,int>, Hang>& H)
{
  cerr << H.first << " -> " << H.second << endl;
  return os;
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::check_mesh2d() const
{
  // check quads

  int cmin = 0, cmax = quads.size(); 
  int vmin = 0, vmax = nnodes(); 

  for(vector<Quad>::const_iterator p = quads.begin();p!=quads.end();p++)
    {
      const Quad& q = *p;

      // check vertex id
      for(int i=0; i< q.nvertexs(); i++)
	{
	  if((q.vertex(i)<vmin)||(q.vertex(i)>vmax))
	    {
	      cerr << "Vertex invalid in Cell: " << *p << " : ";
	      cerr << q.vertex(i) << endl; exit(1);
	    }
	}
      // check child id
      for(int i=0; i< q.nchilds(); i++)
	{
	  if((q.child(i)<cmin)||(q.child(i)>cmax))
	    {
	      cerr << "Chid invalid in Cell: " << *p << " : ";
	      cerr << q.child(i) << endl; exit(1);
	    }
	}
    }

  // check linehang
  LineHangList::const_iterator hp;

  for (hp = LineHang.begin(); hp!=LineHang.end(); ++hp)
  {
    assert(hp->second.rneighbour()>=0);
  }
}

/*---------------------------------------------------*/

int HierarchicalMesh2d::QuadNeighbour(const Quad& q, int e) const
{
  int ie = q.edge(e);
  int m  = edges[ie].master();
  if (q==quads[m]) return edges[ie].slave();
  return m;      
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::init_edges2d()
{
  EdgeManager EM(edges,quads,co2n,eo2n);
  EM.InitEdges();
  EM.SortHangings();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::GetVertexesOfEdge(fixarray<3,int>& v, int e) const
{
  const Edge& E = edge(e);
  const Quad* Q = &quad(E.master());

  int le = E.LocalMasterIndex();
  v[0] = (*Q)[le]; v[1] = (*Q)[(le+1)%4]; v[2] = -1;

  if(!Q->sleep()) 
    {
      if (E.slave()==-1) return;
      Q = &quad(E.slave());
      le = E.LocalSlaveIndex();
    }
  v[2] = QuadLaO.edge_vertex(*Q,le);
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::GetVertexesOfEdge(fixarray<2,int>& v, int e) const
{
  const Edge& E = edge(e);
  const Quad* Q = &quad(E.master());

  int le = E.LocalMasterIndex();
  v[0] = (*Q)[le]; v[1] = (*Q)[(le+1)%4];
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::GetAwakeCells(set<int>& v) const
{
  for (int i=0; i<ncells(); i++)
    {
      if (!quads[i].sleep()) v.insert(i);
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::GetAwakePatchs(set<int>& v) const
{
  for (int i=0; i<ncells(); i++)
    {
      if (!quads[i].sleep()) 
	{
	  int f = quads[i].father();
	  assert(f!=-1);
	  v.insert(f);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::ConstructQ2PatchMesh(IntVector& q2patchmesh) const
{
  typedef IntSet::iterator It;
  IntSet patche;
  GetAwakePatchs(patche);
  vector<IntSet> patch_on_level(nlevels());
  q2patchmesh.resize(0);
  for (It it = patche.begin();it!=patche.end();++it)
    patch_on_level[quads[*it].level()].insert(*it);
  // grobgitterzellen kommen ins patchmesh
  for (It it=patch_on_level[0].begin();it!=patch_on_level[0].end();++it)
    q2patchmesh.push_back(*it);
  // der Rest wird eins groeber
  for (int l=1;l<nlevels();++l)
    {
      It it = patch_on_level[l].begin();
      while (it!=patch_on_level[l].end())
	{
	  int v = quads[*it].father();
	  assert(v!=-1);
	  q2patchmesh.push_back(v);
	  IntVector nk = Nachkommen(v);
	  for (int i=0;i<nk.size();++i)
	    patch_on_level[quads[nk[i]].level()].erase(nk[i]);
	  it = patch_on_level[l].begin();
	}
    }  
}

/*---------------------------------------------------*/

IntVector HierarchicalMesh2d::ConstructQ4Patch(int c) const
{
  IntVector patch(25,-1);

  for(int i=0; i<25; i++)
  {
    // Vertex i steht an Position (x,y)
    int x = i%5;
    int y = i/5;

    // Position von erstem Kind
    int fcx = x/3;
    int fcy = y/3;
    // Index davon
    int fci = fcy*2+abs(fcx-fcy);

    // Position vom Kind im Kind
    int scx = (x-2*fcx)/2;
    int scy = (y-2*fcy)/2;
    // Index davon
    int sci = scy*2+abs(scx-scy);

    // Position des Vertex
    int vx = x-2*fcx-scx;
    int vy = y-2*fcy-scy;
    // Index davon
    int vi = vy*2+abs(vx-vy);

    patch[i] = quads[quads[quads[c].child(fci)].child(sci)].vertex(vi);
  }
  return patch;
}

/*---------------------------------------------------*/


set<int> HierarchicalMesh2d::CellNeighbours(int iq) const
{
  // soll nur wache Zellen zurueckgeben !!
  set<int> neighbors;
  const Quad& q = quad(iq);
  if( q.sleep() ) return neighbors;

  for(int i=0;i<4;i++)
    {
      int ge = q.edge(i);
      const Edge& e = edge(ge);
      int m = e.master();
      int s = e.slave();
      if(m!=-1)
	{
	  if(  !quad(m).sleep() ) neighbors.insert(m); 
	  else
	    {
	      for(int ii=0;ii<4;ii++) neighbors.insert(quad(m).child(ii));
	    }
	}
      // wenn slave groeber, koennen keine Nachbarn gefunden werden !!!
      if(s!=-1)
	{
	  if(  !quad(s).sleep() ) neighbors.insert(s); 
	  else
	    {
	      for(int ii=0;ii<4;ii++) 
		neighbors.insert(quad(s).child(ii));
	    }
	}
    }
  return neighbors;
}

/*------------------------------------------------------*/

int HierarchicalMesh2d::neighbour(int c, int le) const
{
  assert(le<4);
  const Quad& Q = quad(c);
  assert(Q.edge(le)>=0);
  const Edge& E = edge(Q.edge(le));
  int m = E.master();
  int nq = m;
  if (m==c) nq = E.slave();
  return nq;
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::FillAllBoundaryLines()
{
  HangBLList  B;
  // erstmal alle rein
  EdgeVector edge;
  for(int i=0;i<nblines();i++)
    {
      BoundaryLine b = bline(i);
      edge = b.vertex();
      sort(edge.begin(),edge.end());
      B.insert(make_pair(edge,b));
    }
  for(int i=0;i<ncells();i++)
    {
      Quad q = quad(i);
      for(int e=0;e<4;++e)
	{
	  QuadLawOrder().global_edge_unsorted(edge,q,e);
	  sort(edge.begin(),edge.end());
	  HangBLList::iterator p = B.find(edge);
	  if(p==B.end())  
	    {
	      BoundaryLine nb;
	      nb.of_quad()      = i;
	      nb.edge_in_quad() = e;
	      nb.material() = -1;
	      B.insert(make_pair(edge,nb));
	    }
	  else
	    {
	      BoundaryLine& b = p->second;
	      if (b.material()==-1) B.erase(p);
	      else 
		{
		  b.of_quad()      = i;
		  b.edge_in_quad() = e;
		}
	    }
	}
    }
  int n = B.size();
  Blines.resize(n);
  HangBLList::iterator p = B.begin();
  for(int i=0;i<n;i++) 
    {
      Blines[i] = p->second;
      if(Blines[i].material()==-1) 
	{
	  Blines[i].vertex() = p->first;
	}
      p++;
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::WriteAll(const string& name) const
{
  ofstream out(name.c_str());

  out << dimension() << " dimension" << endl;
  out << nnodes() << " vertexs"   << endl;
  out << mnlevels << " mnlevels"   << endl;

  out << "vertexs2d\n" << vertexs2d << endl;
  
  out << "vo2n\n" << vo2n << endl;
  out << "co2n\n" << co2n << endl;
  out << "eo2n\n" << eo2n << endl;
  
  out << "quads\n" << quads << endl;
  out << "Blines\n" << Blines << endl;
  out << "LineHang\n" << LineHang << endl;

  out.close();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::write_inp(const string& bname) const
{
  string name = bname;
  int name_size=name.size();
  if(name_size<4) name += ".inp";
  if(name.substr(name_size-4,4)!=".inp"){
    name += ".inp";
  }

  ofstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh2d::write_inp()\n";
      cerr << "cannot open file " << name << endl;
      abort();
    }
  
  int nt = ncells()+nblines();
  file <<nnodes()<<" "<<nt<<" "<<0<<" "<<0<<" "<<0<<endl;

  
  for(int i=0;i<nnodes();i++) file <<i<<" "<<vertex2d(i) << " " << 0. << endl;

  for(int i=0;i<ncells();i++)
    {
      file << i<<" "<<0<<" quad "<<quad(i).vertex()<<endl;
    }
  for(int i=0;i<nblines();i++)
    {
      file << i<<" "<<bline(i).material()<<" line "<<bline(i).vertex()<<endl;
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::write_vtk(const string& name) const
{
  ofstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh2d::write_vtk()\n";
      cerr << "cannot open file " << name << endl;
      exit(1);
    }

  int nn = nnodes();

  file << "# vtk DataFile Version 2.4 "<<endl;
  file << "output from GascoigneStd" << endl;
  file << "ASCII" << endl;
  file << "DATASET UNSTRUCTURED_GRID" << endl;
  file << "POINTS " << nn << " FLOAT " << endl;

  int ne = ncells();
  for (int i=0; i<nn; i++)
    {
      file<<  vertex2d(i) << " " << 0 << endl;
    }
  
  file << endl << "CELLS " << ne <<" " << 5*ne << endl;
  
  for (int c=0; c<ne; c++)
    {
      file << 4 << " ";
      for(int i=0;i<4;i++)
	{
	  file << vertex_of_cell(c,i) << " "; 
	}
      file << endl; 
    }     
  file << endl << "CELL_TYPES " << ne << endl;
  for(int i=0;i<ne;i++) file << 9 << " ";
}

/*---------------------------------------------------*/

pair<bool,triple<int,int,int> > HierarchicalMesh2d::check_inp(const string& name)
{
  //  detect some errors in input-file... and more

  ifstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh2d::check_inp()\n";
      cerr << "cannot open file " << name << endl;
      exit(1);
    }

  bool first_one = 1;
  int  nv, nl, nq, nt;
  int  n_unkonwn;
  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  Vertex3d c; int ind;
  for(int i=0;i<nv;i++)
    {
      file >> ind >> c;
    }

  nq = 0; nl = 0;
  fixarray<8,int> ih;
  fixarray<4,int> iq;
  fixarray<2,int> il;
  for(int i=0;i<nt;i++)
    {
      string name;
      string mat;
      int ii;
      file >> ii >> mat  >> name;
      if(name=="quad")
	{
	  file >> iq;
	  nq++;
	  if( (iq[0]==0)||(iq[1]==0)||(iq[2]==0)||(iq[3]==0) )
	    {
	      first_one = 0;
	    }
	}
      else if(name=="line")
	{
	  file >> il;
	  nl++;
	  if( (il[0]==0)||(il[1]==0))
	    {
	      first_one = 0;
	    }
	}
    }

  // fehlerabfragen ....
  if(nt!=(nl+nq))
    {
      cerr << "wrong number of cells: " << nt << endl;
      exit(1);
    }

  file.close();

  return make_pair(first_one,make_triple(nl,nq,0));
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::read_inp(const string& name)
{
  // check mesh.....

  pair<bool,tint> p = check_inp(name);
  bool first_one = p.first;
  tint n = p.second;

  int nl = n.first;
  int nq = n.second;

  ifstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh2d::read_inp()\n";
      cerr << "cannot open file " << name << endl;
      abort();
    }

  int  nv,  nt, n_unkonwn;

  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  if( _i_showoutput ) {
    cout << "2D Mesh: " << nv << " nodes, ";
    cout << nl << " lines, " ;
    cout << nq << " quads"<< endl;
  }

  vertexs2d.reserve(nv);
  vertexs2d.resize(nv,Vertex2d());
  quads.    reserve(nq);
  quads.    resize(nq,Quad());
  Blines.   reserve(nl);
  Blines.   resize(nl);

  Vertex2d c;      int ind;  double z;
  for(int i=0;i<nv;i++)
    {
      file >> ind >> c >> z;
      vertexs2d[i] = c;
    }

  fixarray<4,int> iqv;
  fixarray<2,int> ilv;
  int iq = 0;
  int il = 0;
  for(int i=0;i<nt;i++)
    {
      string name;	int unknown; string matstring;
      file >> unknown >> matstring  >> name;
      if(name=="quad")
	{
	  file >> iqv;
	  if(first_one) for(int iii=0;iii<4;iii++) iqv[iii]--;
	  quads[iq++].vertex() = iqv;
	}
      else if(name=="line")
	{
	  file >> ilv;
	  if(first_one) {ilv[0]--; ilv[1]--;}
	  
	  BoundaryLine li;
	  li.material() = atoi(matstring.c_str());
	  li.vertex() = ilv;
	  init_line(li);
	  Blines[il++] = li;
	}
    }
  init_edges2d();
  IntVector nix(0);
  refine(nix,nix);
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::write(const string& bname) const
{
  write_gup(bname);
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::write_gup(const string& bname) const
{
  string name = bname;
  int name_size=name.size();
  if(name_size<4) name += ".gup";
  if(name.substr(name_size-4,4)!=".gup"){
    name += ".gup";
  }

  ofstream out(name.c_str());

  out.precision(16);
  out << dimension() << " dimension" << endl;
  out << nnodes() << " vertexs"   << endl;

  for (int i=0; i<nnodes(); i++)
    {
      out << " " << vertex2d(i) << endl;
    }
  out << quads.size() << " quads" << endl;
  
  for (int i=0; i<quads.size(); i++)
    {
      out << quad(i);
    }
  out << LineHang << endl;
  
  out << Blines.size() << " boundarylines" << endl;
  for (int i=0; i<Blines.size(); i++)
    {
      out << Blines[i].material() << " " << Blines[i];
    }
  out << endl << endl << edges.size() << " edges" << endl;
  for (int i=0; i<edges.size(); i++)
    {
      out << " " << edges[i];
    }
  out << endl;
  out.close();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::write_gip(const string& bname) const
{
  string name = bname;
  int name_size=name.size();
  if(name_size<4) name += ".gip";
  if(name.substr(name_size-4,4)!=".gip"){
    name += ".gip";
  }

  ofstream out(name.c_str(), ios_base::out|ios_base::binary);

  out.precision(16);
  int dim=dimension(), n=nnodes(), sizeInt=sizeof(int);
  out.write(reinterpret_cast<const char*>(&dim),sizeInt);
  out.write(reinterpret_cast<const char*>(&n),sizeInt);

  for (int i=0; i<nnodes(); i++)
    {
      vertex2d(i).BinWrite(out);
    }

  int nquads = quads.size();
  out.write(reinterpret_cast<const char*>(&nquads),sizeInt);
 
  for (int i=0; i<quads.size(); i++)
    {
      quad(i).BinWrite(out);
    }

  LineHang.BinWrite(out);
  int nblines=Blines.size();
  out.write(reinterpret_cast<const char*>(&nblines),sizeInt);
  for (int i=0; i<Blines.size(); i++)
    {
      int mat = Blines[i].material();
      out.write(reinterpret_cast<const char*>(&mat),sizeInt);
      Blines[i].BinWrite(out);
    }
  int nedges=edges.size();
  out.write(reinterpret_cast<const char*>(&nedges),sizeInt);
  for (int i=0; i<edges.size(); i++)
    {
      edges[i].BinWrite(out);
    }
  out.close();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::read_gup(const string& bname)
{
  string name = bname;
  int name_size=name.size();
  if(name_size<4) name += ".gup";
  if(name.substr(name_size-4,4)!=".gup"){
    name += ".gup";
  }

  string symbol;

  ifstream file(name.c_str());

  if(!file.is_open())
    {
      cerr << "HierarchicalMesh2d::read_gup(): error in file "<< name<<endl;
      abort();
    }

  Blines.clear();
  edges.clear();
  LineHang.clear();

  int n, dim;
  file >> dim >> symbol >> n >> symbol;

  assert(dim==2);

  if( _i_showoutput ){
    cout << "Mesh 2d  : ";
  }
  vertexs2d.reserve(n);
  vertexs2d.resize(n);

  if (symbol!="vertexs")
    {
      cout << "HierarchicalMesh2d::read error" << endl;
      exit(1);
    }
  for (int i=0; i<n; i++)
    {
      file >> vertexs2d[i];
    }
  if( _i_showoutput ){
    cout << n << " nodes, ";
  }

  file >> n >> symbol;

  if (symbol!="quads")
    {
      cout << "HierarchicalMesh2d::read error 2" << endl;
      exit(1);
	}
  quads.reserve(n);
  quads.resize(n);
  for (int i=0; i<quads.size(); i++)
    {
      file >> quads[i];
    }
  if( _i_showoutput ){
    cout <<  n << " quads, ";
  }
  file >> LineHang;
  if( _i_showoutput ){
    cout << LineHang.size() << " hangs, ";
  }
  
  file >> n >> symbol;
  int number = 0;
  if (symbol=="boundary_materials")
    {
      /* old format */
      BoundaryLine bol;
      for (int i=0; i<n; i++)
	{
	  file >> symbol;
	  if (symbol!="material")
	    {
	      cout << "HierarchicalMesh2d::read error 4" << endl;
	      exit(1);
	    }
	  int nn;
	  file >> bol.material() >> nn;
	  for (int j=0; j<nn; j++)
	    {
	      file >> bol;
	      Blines.push_back(bol);
	    } 
	  number += nn;
	}
    }
  else if (symbol=="boundarylines")
    {
      /* new format */
      BoundaryLine bol;
      for (int i=0; i<n; i++)
	{
	  file >> bol.material() >> bol;
	  Blines.push_back(bol);
	}
      number = n;
    }
  else
    {
      cout << "HierarchicalMesh2d::read error 3" << endl;
      exit(1);
    }
  if( _i_showoutput ){
    cout << number << " lines, ";
  }
  file >> n >> symbol;
  if( _i_showoutput ){
    cout << n << " edges" << endl;
  }
  if (symbol=="edges")
    {
      Edge e;
      for (int i=0; i<n; i++)
	{
	  file >> e;
	  edges.push_back(e);
	}
    }
  if (edges.size()==0) init_edges2d();
  post_refine2d();
  file.close();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::read_gip (const string& bname)
{
  string name = bname;
  int name_size=name.size();
  if(name_size<4) name += ".gip";
  if(name.substr(name_size-4,4)!=".gip"){
    name += ".gip";
  }

  string symbol;

  ifstream file(name.c_str(), ios_base::in|ios_base::binary);

  if(!file.is_open())
    {
      cerr << "HierarchicalMesh2d::read_gip(): error in file "<< name <<endl;
      abort();
    }

  Blines.clear();
  edges.clear();
  LineHang.clear();

  int n, dim, sizeInt=sizeof(int);
  file.read(reinterpret_cast<char*>(&dim),sizeInt);
  file.read(reinterpret_cast<char*>(&n),sizeInt);

  assert(dim==2);

  if( _i_showoutput ){
    cout << "Mesh 2d  : ";
  }
  vertexs2d.reserve(n);
  vertexs2d.resize(n);

  for (int i=0; i<n; i++)
    {
      vertexs2d[i].BinRead(file);
    }
  if( _i_showoutput ){
    cout << n << " nodes, ";
  }
  file.read(reinterpret_cast<char*>(&n),sizeInt);
  quads.reserve(n);
  quads.resize(n);
  for (int i=0; i<quads.size(); i++)
    {
      quads[i].BinRead(file);
    }
  if( _i_showoutput){
    cout <<  n << " quads, ";
  }
  LineHang.BinRead(file);
  if( _i_showoutput){
    cout << LineHang.size() << " hangs, ";
  }
  file.read(reinterpret_cast<char*>(&n),sizeInt);
  int number = 0;
  BoundaryLine bol;
  for (int i=0; i<n; i++)
    {
      file.read(reinterpret_cast<char*>(&bol.material()),sizeInt);
      bol.BinRead(file);
      Blines.push_back(bol);
    }
  number = n;
  if( _i_showoutput){
    cout << number << " lines, ";
  }
  file.read(reinterpret_cast<char*>(&n),sizeInt);
  if( _i_showoutput){
    cout << n << " edges" << endl;
  }
  Edge e;
  for (int i=0; i<n; i++)
    {
      e.BinRead(file);
      edges.push_back(e);
    }
  if (edges.size()==0) init_edges2d();
  post_refine2d();
  file.close();
}

/*---------------------------------------------------*/

int HierarchicalMesh2d::regular_grid2d_one(IntSet& celllist, 
					 IntVector& coarsesub,
					 IntSet& CellRefList,
					 IntSet& CellCoarseList) const
{
  /* detects jump over two levels across LineHangs */

  int n = 0;
  LineHangList::const_iterator  hp = LineHang.begin();
 
  for (; hp!=LineHang.end(); ++hp)
  {
    int cr = hp->second.rneighbour();
    int cn = hp->second.cneighbour();

    assert(cr>=0);

    if(cn!=-1)
      {
	if(quad(cr).childs().size()==0) continue;

	fixarray<2,int>  f;
	QuadLaO.childs_of_global_edge(f,quad(cr),hp->first);

	for(unsigned i=0;i<f.size();++i)
	  {
	    int c = f[i];
	    if ( (CellRefList.find(c)!=CellRefList.end()) || 
		 (quads[c].sleep()) )
	      {
		if (CellCoarseList.find(cn)!=CellCoarseList.end())
		  {
		    coarsesub.push_back(cn);
		    break;
		  }
		else
		  {
		    assert(!quad(cn).sleep());

		    pair<IntSetIt,bool> p = celllist.insert(cn);
		    if(p.second)  { n++; break; }
		  }
	      }
	  }
      }
  }
  //unique(coarsesub.begin(),coarsesub.end());
  return n+coarsesub.size();
}

/*---------------------------------------------------*/

int HierarchicalMesh2d::regular_grid2d_two(IntSet& celllist,
				       IntSet& CellRefList) const
{
  /* detects more than 2 HangLines on one Quad */

  int nto = 0;
  IntVector nh(quads.size());
  LineHangList::const_iterator p = LineHang.begin();
  for (; p!=LineHang.end(); p++)
    {
      int i = p->second.cneighbour();
      if (i<0)                                    continue;
      if (quad(i).sleep())                        continue;
      if (CellRefList.find(i)!=CellRefList.end()) continue;

      nh[i]++;
    }
  for (int i=0; i<quads.size(); i++)
    {
      if (nh[i]>2)
	{
	  pair<IntSetIt,bool> pp = celllist.insert(i);
	  if(pp.second)  nto++;
	}
    }
  return nto;
}

/*---------------------------------------------------*/

int HierarchicalMesh2d::smooth_edges()
{
  IntSet         CellRefList,CellCoarseList,coarseadd;
  HangContainer2d hangset(LineHang);

  int r = 0;//regular_grid2d_three(CellRefList);

  ghost2d(hangset,CellRefList,coarseadd);

  basic_refine2d(hangset,CellRefList,CellCoarseList);
  post_refine2d();

  return r;
}

/*---------------------------------------------------*/

int HierarchicalMesh2d::regular_grid2d_three(IntSet& CellRef,
					     IntSet& CellCoarse) const
{
  vector<int> maxlevel(nnodes());
  vector<int> minlevel(nnodes(),1000);

  // set maximal levels for vertices
  for (int i=0; i<quads.size(); i++)
    {
      const Quad& q = quad(i);
      if(q.sleep()) continue;

      int lev = q.level();
      if (CellRef.find(i)!=CellRef.end()) lev++;
      
      int father = q.father();
      if (father>=0)
	{
	  if (CellCoarse.find(father)!=CellCoarse.end())
	    {
	      lev--;
	    }
	}
      for (int j=0; j<4; j++)
	{
	  int k = q[j];
	  maxlevel[k] = Gascoigne::max_int( maxlevel[k], lev);
	  minlevel[k] = Gascoigne::min_int( minlevel[k], lev);
	}
    }
  set<int> cand;
  for (int i=0; i<nnodes(); i++)
    {
      if (maxlevel[i]>=minlevel[i]+2)
	{
	  cand.insert(i);
	}
    }

  int ref = 0;
  for (int i=0; i<quads.size(); i++)
    {
      const Quad& q = quad(i);
      if (q.sleep()) continue;
      int lev = q.level();
      for (int j=0; j<4; j++)
	{
	  int k = q[j];
	  set<int>::const_iterator p = cand.find(k);
	  if (p!=cand.end())
	    {
	      if (CellRef.find(i)==CellRef.end())
		{
		  if (maxlevel[k]>=lev+2) CellRef.insert(i);
		  ref++;
		}
	    }
	}
    }
  IntSet coarsesub;
  IntSet::const_iterator p = CellCoarse.begin();
  for ( ; p!= CellCoarse.end(); p++)
    {
      const Quad& q = quad(*p);
      int lev = q.level();
      for (int j=0; j<4; j++)
	{
	  int k = q[j];
	  if (lev<=minlevel[k]-1)
	    {
	      coarsesub.insert(*p);
	    }
	}
    }
  int coarse = 0;
  p = coarsesub.begin();
  for ( ; p!= coarsesub.end(); p++)
    {
      if (CellCoarse.find(*p)!=CellCoarse.end())
	{
	  CellCoarse.erase(*p);
	  coarse++;
	}
    }
  //cout << "(" << ref << " " << coarse << ")\n";
  return ref+coarse;
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::GetMinMaxLevels(IntVector& maxi, IntVector& mini,
					 const IntSet& CellRef) const
{
  // set maximal levels for vertices
  //
  maxi.resize(nnodes());
  mini.resize(nnodes()); mini = 1000;
  for (int i=0; i<quads.size(); i++)
    {
      const Quad& q = quad(i);
      if(q.sleep()) continue;

      int lev = q.level();
      if (CellRef.find(i)!=CellRef.end()) lev++;
      for (int j=0; j<4; j++)
	{
	  int k = q[j];
	  maxi[k] = Gascoigne::max_int( maxi[k], lev);
	  mini[k] = Gascoigne::min_int( mini[k], lev);
	}
    }
}

/*---------------------------------------------------*/

int HierarchicalMesh2d::regular_grid2d_three_refine(IntSet& CellRef) const
{
  IntVector maxlevel, minlevel;

  GetMinMaxLevels(maxlevel,minlevel,CellRef);
  
  set<int> cand;
  for (int i=0; i<nnodes(); i++)
    {
      if (maxlevel[i]>=minlevel[i]+2)
	{
	  cand.insert(i);
	}
    }

  int ref = 0;
  for (int i=0; i<quads.size(); i++)
    {
      const Quad& q = quad(i);
      if (q.sleep()) continue;
      int lev = q.level();
      for (int j=0; j<4; j++)
	{
	  int k = q[j];
	  set<int>::const_iterator p = cand.find(k);
	  if (p!=cand.end())
	    {
	      if (CellRef.find(i)==CellRef.end())
		{
		  if (maxlevel[k]>=lev+2) CellRef.insert(i);
		  ref++;
		}
	    }
	}
    }
  //cout << "(" << ref << ")";
  return ref;
}

/*---------------------------------------------------*/

int HierarchicalMesh2d::regular_grid2d_three_coarse(IntSet& CellRef,
						    IntSet& CellCoarse) const
{
  int maxl = 0;
  vector<IntSet> LevelCellCoarse;
  {
    IntSet::const_iterator p = CellCoarse.begin();
    for ( ; p!= CellCoarse.end(); p++)
      {
	const Quad& q = quad(*p);
	int lev = q.level();
	maxl = Gascoigne::max_int(maxl,lev);
      }
    LevelCellCoarse.resize(maxl+1);
    p = CellCoarse.begin();
    for ( ; p!= CellCoarse.end(); p++)
      {
	const Quad& q = quad(*p);
	int lev = q.level();
	LevelCellCoarse[lev].insert(*p);
      }
  }
  int coarse = 0;
  for (int i=maxl; i>=0; i--)
    {
      IntVector maxlevel, minlevel;

      GetMinMaxLevels(maxlevel,minlevel,CellRef);
  
      IntSet coarsesub;
      IntSet::const_iterator p = LevelCellCoarse[i].begin();
      for ( ; p!= LevelCellCoarse[i].end(); p++)
	{
	  const Quad& q = quad(*p);
	  int lev = q.level();
	  for (int j=0; j<4; j++)
	    {
	      int k = q[j];
	      if (lev+1<maxlevel[k])
		{
		  coarsesub.insert(*p);
		  continue;
		}
	    }
	}
      //CellCoarse.erase(coarsesub.begin(),coarsesub.end());
      p = coarsesub.begin();
      for ( ; p!= coarsesub.end(); p++)
	{
	  if (CellCoarse.find(*p)!=CellCoarse.end())
	    {
	      CellCoarse.erase(*p);
	      coarse++;
	    }
	}
    }
  //cout << "(" << coarse << ")";
  return coarse;
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::global_coarse()
{
  IntSet RefList,CoarseList;

  for (int i=0; i<quads.size(); i++)
    {
      const Quad& q = quads[i];
      if ( q.sleep()) continue;
      if (!q.level()) continue;
      CoarseList.insert(q.father());
    }

  HangContainer2d hangset(LineHang);

  ghostglobalcoarse(hangset,CoarseList);

  basic_refine2d(hangset,RefList,CoarseList);

  post_refine2d();
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::LoadFathers(IntVector& v) const
{
  IntSet fathers;
  for (int i=0; i<v.size(); i++)
    {
      if (quad(v[i]).level()==0) continue;
      int f = quad(v[i]).father();
      if (f<0) 
	{
	  cerr << "HierarchicalMesh2d::LoadFathers no father\n";
	  abort();
	}
      fathers.insert(f);
    }
  Set2Vec(v,fathers);
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::FillVertexLevels(IntVector& dst) const
{
  dst.resize(nnodes(),100);
  for (int i=0; i<ncells(); i++)
    {
      const Quad& Q = quads[i];
      if (Q.sleep()) continue;

      int level = Q.level();
      for (int j=0; j<4; j++)
	{
	  int k = Q[j];
	  dst[k] = Gascoigne::min_int(dst[k],level);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::RefineCoarseNodes
(IntSet& dst, const IntVector& refnodes, const IntVector& vertexlevel) const
{
  IntSet h;

  Vec2Set(h,refnodes);
  
  dst.clear();
  for (int i=0; i<ncells(); i++)
    {
      const Quad& Q = quads[i];
      if (Q.sleep()) continue;

      int f = Q.father();
      if (f<0)  continue;

      const Quad& QF = quads[f];
	   
      for (int j=0; j<4; j++)
	{
	  int k = Q[j];
	  if (h.find(k)==h.end()) continue;

	  int minlevel = vertexlevel[QF[0]];
	  for (int v=1; v<4; v++)
	    {
	      minlevel = Gascoigne::min_int(minlevel,vertexlevel[QF[v]]);
	    }
	  for (int v=0; v<4; v++)
	    {
	      int w = QF[v];
	      assert(w<vertexlevel.size());
	      if (vertexlevel[w]==minlevel)
		{
                  dst.insert(w);
		}
	    }
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::VertexToCells(IntVector& dst, const IntSet& src, 
				       const IntVector& vertexlevel) const
{
  for (int i=0; i<ncells(); i++)
    {
      const Quad& Q = quads[i];
      if (Q.sleep()) continue;
      int level = Q.level();
      for (int j=0; j<4; j++)
	{
	  int k = Q[j];
	  if (vertexlevel[k]==level)
	    {
	      if (src.find(k)!=src.end())
		{
		  dst.push_back(i);
		  break;
		}
	    }
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::VertexToCellsCoarsening(IntVector& dst, const IntSet& src, 
						 const IntVector& vertexlevel) const
{
  int limit = 2;

  for (int i=0; i<ncells(); i++)
    {
      const Quad& Q = quads[i];
      if (Q.sleep()) continue;
      int count = 0;
      for (int j=0; j<4; j++)
	{
	  int k = Q[j];
	  if (src.find(k)!=src.end())
	    {
	      count++;
	    }
	}
      if (count>=limit) dst.push_back(i);
    }
}

/*---------------------------------------------------*/

// void HierarchicalMesh2d::double_patch_refine
// (IntVector& cell_ref, IntVector& cell_coarse)
// {
//   LoadFathers(cell_ref);
//   LoadFathers(cell_ref);
//   LoadFathers(cell_coarse);
//   LoadFathers(cell_coarse);

//   DoubleCoarseHierarchicalMesh2d CM(*this);

//   CM.BasicInit();
//   CM.refine(cell_ref,cell_coarse);
//   CM.GetRefinedList(cell_ref);
//   CM.GetCoarsedList(cell_coarse);

//   IntVector ref(0), coarse(0);

//   for (int i=0; i<cell_ref.size(); i++)
//     {
//       const Quad& q = quad(cell_ref[i]);
//       if (!q.sleep())
// 	{
// 	  cout << "HierarchicalMesh2d::patchrefine am dampfen1 !";
// 	  abort();
// 	}
//       for (int j=0; j<4; j++)
// 	{
// 	  int child = q.child(j);
// 	  for (int jj=0; jj<4; jj++)
// 	    {
// 	      int grandchild = quad(child).child(jj);
// 	      ref.push_back(grandchild);
// 	    }
// 	}
//     }  
//   cell_coarse.resize(0);
//   refine(ref,cell_coarse);
// }

/*---------------------------------------------------*/

void HierarchicalMesh2d::recursive_childs(int q, IntVector& ref, int d) const
{
  if (d>0)
    {
      const Quad& Q = quad(q);
      assert(Q.sleep());
      for (int j=0; j<4; j++)
	{
	  int child = Q.child(j);
	  recursive_childs(child,ref,d-1);
	}
    }
  else
    {
      ref.push_back(q);      
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::patch_refine
(IntVector& cell_ref, IntVector& cell_coarse)
{
  for (int i=0; i<pdepth; i++)
    {
      LoadFathers(cell_ref);
      LoadFathers(cell_coarse);
    }

  CoarseHierarchicalMesh2d CM(*this);

  CM.BasicInit(pdepth);
  CM.refine(cell_ref,cell_coarse);
  CM.GetRefinedList(cell_ref);
  CM.GetCoarsedList(cell_coarse);

  IntVector ref(0), coarse(0);
 
  for (int i=0; i<cell_ref.size(); i++)
    {
      recursive_childs(cell_ref[i],ref,pdepth);
    }  
  for (int i=0; i<cell_coarse.size(); i++)
    {
      recursive_childs(cell_coarse[i],coarse,pdepth+1);
    } 
  refine(ref,coarse);
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::FillVolumes(DoubleVector& vol) const
{
  vol.resize(ncells(),0.);
  for (int i=0; i<ncells(); i++)
    {
      const Quad& Q = quad(i);
      Vertex2d V = vertex2d(Q[2]);
      V -= vertex2d(Q[0]);
      vol[i] = V*V;
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh2d::writeq2(const IntVector &a,const vector<int> & b,int np) const
{
  char s[29];
  for (int p=0;p<np;++p)
    {
      sprintf(s,"n_%d_%d",ncells(),p);
      ofstream aus(s);

      for (int i=0; i<a.size(); ++i)
	{
	  const Quad& Q = quads[a[i]];
	  if (b[i]==p)
	    {
	      aus << vertex2d(Q[0]) << endl;
	      aus << vertex2d(Q[1]) << endl;
	      aus << vertex2d(Q[2]) << endl;
	      aus << vertex2d(Q[3]) << endl;
	      aus << vertex2d(Q[0]) << endl << endl;
	    }
	}
      aus.close();
    }
}
}

/*---------------------------------------------------*/

