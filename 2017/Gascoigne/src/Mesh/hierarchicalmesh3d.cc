/**
*
* Copyright (C) 2004, 2005, 2006, 2008, 2010, 2011 by the Gascoigne 3D authors
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


#include  "hierarchicalmesh3d.h"
#include  "deletecells.h"
#include  "vecalgo.h"
#include  "set2vec.h"
#include  "quad.h"
#include  "levelcomparer3d.h"
#include  "facemanager.h"
#include  "stlio.h"
#include  "regular_update.h"
#include  "coarsehierarchicalmesh3d.h"

#include  <fstream>
#include  "giota.h"


using namespace std;

namespace Gascoigne
{
typedef  triple<int,int,int>          tint;

/*------------------------------------------------------*/

HierarchicalMesh3d::HierarchicalMesh3d() 
  : HierarchicalMesh(), HexLaO(hexs) {}

/*------------------------------------------------------*/

HierarchicalMesh3d::HierarchicalMesh3d(const ParamFile* paramfile) 
  : HierarchicalMesh(), HexLaO(hexs)
{
  BasicInit(paramfile);
}

/*------------------------------------------------------*/

HierarchicalMesh3d::HierarchicalMesh3d(const HierarchicalMesh3d& H)
  : HierarchicalMesh(), HexLaO(hexs)
{
  *this = H;
}

/*------------------------------------------------------*/

HierarchicalMesh3d& HierarchicalMesh3d::operator=(const HierarchicalMesh3d& H)
{
  HierarchicalMesh::operator= (H);
  // copy all data
  vertexs3d = H.vertex3d();
  hexs      = H.hex();
  Bquads    = H.bquad();
  LineHang  = H.linehang();
  QuadHang  = H.quadhang();
  hexofcurved = H.GetHexOfCurved();

  return *this;
}

/*------------------------------------------------------*/

pair<int,int> HierarchicalMesh3d::GetBoundaryInformation(int i) const
{
  int material = -1; 
  int le = -1;
  int ib = GetBoundaryCellOfCurved(i);
  if (ib>=0)
    {
      material = bquad(ib).material();
      le       = bquad(ib).edge_in_quad();
    }
  return make_pair(material,le);
}

/*------------------------------------------------------*/

const BoundaryFunction<3>* HierarchicalMesh3d::quad_shape(int i) const
{ 
  if (GetCurvedShapes().empty()) return NULL;

  if (GetCurvedShapes().Curved(i)) return &(GetCurvedShapes().GetShape(i));
 
  return 0;
}

/*------------------------------------------------------*/

set<int> HierarchicalMesh3d::GetColors() const
{
  set<int> coleur;

  for(int i=0;i<nbquads();i++)
    {
      coleur.insert(bquad(i).material());
    }
  return coleur;
}

/*------------------------------------------------------*/

int HierarchicalMesh3d::FindPatchDepth() const
{
  //simple version, sucht nur p=1, p=0
  for(int i=0;i<ncells();i++)
    {
      const Hex& q = hex(i);
      if(q.sleep()) continue;
      int father = q.father();
      if(father==-1) return 0;
      const Hex& qf = hex(father);
      for(int ii=0;ii<8;ii++) 
	{
	  if(hex(qf.child(ii)).sleep()) return 0;
	}
    }
  return 1;
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::InitHexOfCurved()
{
  hexofcurved.clear();

  if (GetCurvedShapes().empty()) return;

  for(int il=0;il<nbquads();++il)
    {
      const BoundaryQuad& B = bquad(il);
      if (GetCurvedShapes().Curved(B.material()))
	{
	  int iq = B.of_quad();
	  hexofcurved.insert(make_pair(iq,il));
	}
    }
}

/*------------------------------------------------------*/

void HierarchicalMesh3d::prepare3d(const IntVector& cell_ref, 
				 const IntVector& cell_coarse,
				 IntSet& CellRefList, 
				 IntSet& CellCoarseList)
{
  /* copies cell_ref into CellRefList without duplets  */

  for (IntVector::const_iterator cp=cell_ref.begin(); 
       cp!=cell_ref.end(); ++cp)
    {
      int c = *cp;
      if ( (c>=0) && (c<hexs.size()) )
	{
	  if(!hexs[c].sleep())
	    {
	      CellRefList.insert(c);
	    }
	}
    }

  /* copies cell_coarse into CellCoarseList without duplets  
     checks if coarse cell in refine list */

  IntSet  help;  

  for (int i=0; i<cell_coarse.size(); i++) 
    {
      int ic = cell_coarse[i];
      if((ic<0)||(ic>=hexs.size())) continue;

      if (  hexs[ic].sleep() )                       continue;
      if ( !hexs[ic].level() )                       continue;
      if ( CellRefList.find(ic)!=CellRefList.end() ) continue;

      help.insert(ic);
    }

  /* checks if coarse cell is cneighbour of a hang */

  for(HangList<4>::const_iterator Lp = QuadHang.begin();
      Lp!=QuadHang.end();Lp++)
    {
      int cn = Lp->second.cneighbour();
      if( help.find(cn)!=help.end() )    help.erase(cn);
    }

  /* marks the father */

  multiset<int> ff;
  
  for (IntSet::const_iterator hp = help.begin();
       hp != help.end(); ++hp)
    {
      ff.insert(hex(*hp).father());
    }

  for (multiset<int>::iterator fp = ff.begin();
       fp != ff.end(); ++fp)
    {
      if (ff.count(*fp)==8) 
	{
	  CellCoarseList.insert(*fp);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::FaceCoarse(HangContainer3d& hangset,
				    const IntSet& cellcoarse) const
{
  fixarray<4,int> quadglob;
  for (IntSetIt cp=cellcoarse.begin(); cp!=cellcoarse.end(); ++cp)
    {
      int f = *cp;
      for(int face=0; face<6; ++face)
	{
	  HexLaO.global_face_unsorted(quadglob,hex(f),face);

	  int ve = HexLaO.face_vertex(hex(f),face);
	  hangset.face_coarse(quadglob,f,ve);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::FaceRefine(HangContainer3d& hangset,
				    const IntSet& cellref) const
{
  fixarray<4,int> quadglob;
  for (IntSetIt cp=cellref.begin(); cp!=cellref.end(); ++cp)
    {
      int f = *cp;
      for(int face=0; face<6; ++face)
	{
	  HexLaO.global_face_unsorted(quadglob,hex(f),face);
	  hangset.face_refine(quadglob,f);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::UpdateHangingEdges(HangContainer3d& hangset,
					    const IntSet& cellref, 
					    const IntSet& cellcoarse) const
{
  HangList<2> oldhangs(LineHang);
  hangset.build_hanging_lines(oldhangs);
  fixarray<2,int> lineglob;
  
  for (IntSetIt cp=cellcoarse.begin(); cp!=cellcoarse.end(); ++cp)
    {
      int f = *cp;
      for(int edge=0; edge<12; ++edge)
	{
	  HexLaO.global_edge_unsorted(lineglob,hex(f),edge);

	  int ve = HexLaO.edge_vertex(hex(f),edge);
	  hangset.line_coarse(lineglob,f,ve);
	}
    }
  for (IntSetIt cp=cellref.begin(); cp!=cellref.end(); ++cp)
    {
      int f = *cp;
      for(int edge=0; edge<12; ++edge)
	{
	  HexLaO.global_edge_unsorted(lineglob,hex(f),edge);

	  hangset.line_refine(lineglob,f,oldhangs);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::UpdateHangs(HangContainer3d& hangset,
				     const IntSet& cellref, 
				     const IntSet& cellcoarse)
{
  // faces
  FaceCoarse(hangset,cellcoarse);
  FaceRefine(hangset,cellref);
  ghost_fill_neighbours3d();

  // edges
  UpdateHangingEdges(hangset,cellref,cellcoarse);
  //  ghost_fill_neighbours2d();
  hangset.make_consistent();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::ghost_fill_neighbours3d()
{
  FaceVector quadglob;
  for(int ic=0;ic<hexs.size();ic++)
    {
      for(int i=0;i<6;i++)
	{
	  HexLaO.global_face_unsorted(quadglob,hex(ic),i);
	  HangList<4>::iterator p = QuadHang.find(quadglob);

	  if(p!=QuadHang.end())
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

void HierarchicalMesh3d::ghost_fill_neighbours2d()
{
  EdgeVector lineglob;
  for(int ic=0; ic<hexs.size(); ic++)
    {
      for(int i=0; i<12; i++)
	{
	  HexLaO.global_edge_unsorted(lineglob,hex(ic),i);
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

void HierarchicalMesh3d::basic_refine3d(HangContainer3d& hangset,
				      const IntSet& CellRefList,
				      const IntSet& CellCoarseList)
{
  int ov = nnodes();
  int oc = hexs.size();
  int oe = edges.size();

  int csub = 8*CellCoarseList.size();
  int cadd = 8*CellRefList.size();

  int vsub = hangset.nDel()+CellCoarseList.size();
  int vadd = hangset.nNew()+CellRefList.size();

  int cdiff = cadd - csub;
  int vdiff = vadd - vsub;

  int nv = ov + vdiff;
  int nc = oc + cdiff;

  clear_transfer_lists();

  /* delete vertex */

  IntVector  vdel, cdel,edel;

  for(IntSetCIt p=CellCoarseList.begin();
      p!=CellCoarseList.end();p++)
    {
      const Hex& q = hex(*p);
      vdel.push_back(HexLaO.middle_vertex(q));
      cdel.insert(cdel.end(),q.childs().begin(),q.childs().end());
    }
  hangset.load_elimination(vdel);
  hangset.NeighbourSwapper();

  //////
  FaceManager EM(edges,hexs,co2n,eo2n);
  //////

  if (withfaces) EM.LoadFaceElimination(edel,CellCoarseList,hangset);
 
  IntSet  QuadRefList, QuadCoarseList, ccdel;

  boundary_prepare3d(QuadRefList, QuadCoarseList, ccdel, hangset);

  transfer(oc,co2n,cdel);
  transfer(ov,vo2n,vdel);
  transfer(oe,eo2n,edel);

  IntVector cnew(cadd), vnew(vadd);
  iota(cnew.begin(),cnew.end(),oc-csub);
  iota(vnew.begin(),vnew.end(),ov-vsub);

  delete_vertexs3d(vo2n);

  if (withfaces) EM.DeleteFaces();

  delete_cells<Hex>(CellCoarseList,hexs,co2n,vo2n);

  vertexs3d.reserve(nv); vertexs3d.resize(nv);
  hexs     .reserve(nc); hexs     .resize(nc);

  hangset.update_olds(vo2n,co2n);
  hangset.update_news(vnew,CellRefList.size());

  new_vertexs3d(hangset,vnew,CellRefList);
  new_hexs     (hangset,cnew,vnew,ov,CellRefList);

  if (withfaces) EM.Check(hangset);

  basic_fill_neighbours3d();

  new_boundary3d(QuadRefList, QuadCoarseList, ccdel);

  IntSet adjustvertex;
  boundary_newton3d(adjustvertex);
  inner_vertex_newton3d(vnew,CellRefList,adjustvertex);

  if (withfaces) EM.Build(CellRefList,hangset);
  Testing();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::Testing()
{
  for(int i=0;i<edges.size();i++)
    {
      int m = edges[i].master();
      int s = edges[i].slave();

      if(s<0) continue;
      if (!hex(m).sleep() && hex(s).sleep())
	{ 
	  swap(edges[i].master(),edges[i].slave());
	  swap(edges[i].LocalMasterIndex(),edges[i].LocalSlaveIndex());
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::boundary_newton3d(IntSet& adjustvertex)
{
  if (GetCurvedShapes().empty()) return;

  for (int i=0; i<Bquads.size(); i++)
    {
      BoundaryQuad& bl = Bquads[i];
      int color = bl.material();
      
      if (GetCurvedShapes().Curved(color))
	{
	  for (int j=0; j<4; j++)
	    {
	      int k = bl.vertex(j);

	      GetCurvedShapes().newton(color, vertexs3d[k]);
	      adjustvertex.insert(k);
	    }
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::boundary_prepare3d(IntSet& QuadRefList, 
					  IntSet& QuadCoarseList,
					  IntSet& ccdel,
					  const HangContainer3d& hangset)
{
  for (int i=0; i<Bquads.size(); i++)
    {
      const BoundaryQuad& bl = Bquads[i];
      
      FaceVector   lineglob;
      lineglob[0] = bl.vertex(0);
      lineglob[1] = bl.vertex(1);
      lineglob[2] = bl.vertex(2);
      lineglob[3] = bl.vertex(3);
      sort(lineglob.begin(),lineglob.end());
      
      if (bl.sleep())
	{
	  if (hangset.ToBeDeleted(lineglob))
	    {
	      QuadCoarseList.insert(i);
	      for (int j=0; j<4; j++) ccdel.insert(bl.child(j));
	    }
	}
      else
	{
	  if (hangset.ToBeCreated(lineglob))
	    {
	      QuadRefList.insert(i);
	    }
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::new_boundary3d(IntSet& ref, IntSet& coarse,
				      IntSet& ccdel)
{
  int oc   = Bquads.size();
  int csub = 4 * coarse.size();
  int cadd = 4 * ref.size();
  int nc   = oc+cadd-csub;
  
  IntVector lo2n;
  transfer(oc,lo2n,ccdel);
  delete_cells<BoundaryQuad>(coarse,Bquads,lo2n,vo2n);
  
  update_boundary_data3d(coarse);
  
  IntVector cnew(cadd);
  iota(cnew.begin(),cnew.end(),oc-csub);
  
  Bquads.reserve(nc);
  Bquads.resize (nc);
  
  new_bquads(lo2n,cnew,ref);
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::update_boundary_data3d(const IntSet& LCoarse)
{
  int no = Bquads.size() - 4*LCoarse.size();
  for(int i=0; i<no; ++i)
    {
      int oq = Bquads[i].of_quad();

      assert(co2n[oq]>=0);

      Bquads[i].of_quad() = co2n[oq];
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::new_bquads(const IntVector& lo2n, 
				  const IntVector& cnew, const IntSet& LRef) 
{
  int nci = 0;

  vector<Quad> emptyq;

  for(IntSetIt cp=LRef.begin(); cp!=LRef.end(); ++cp)
    {
      int father = lo2n[*cp];

      assert(father!=-1);

      BoundaryQuad& bqf = Bquads[father];
      // change father boundary
      vector<int>& qc = bqf.childs();
      qc.resize(4);

      int face = bqf.edge_in_quad();  // face index in hex
      int hexi = bqf.of_quad();       // hex index of father
      FaceVector chvec;
      HexLaO.childs_of_face(chvec,hex(hexi),face);

      for(int ic=0;ic<4;ic++)
	{
	  int childhex = chvec[ic];
	  int inold = cnew[nci+ic];
	  // set childs in father
	  qc[ic]    = inold;
	  // set properties of childs
	  BoundaryQuad& bl = Bquads[inold];
	  bl.level   () = bqf.level()+1;
	  bl.father  () = father;
	  bl.material() = bqf.material();
	  bl.childs().resize(0);

	  HexLaO.load_face(bl.vertex(),hex(childhex),face);
	  bl.of_quad()      = childhex;
	  bl.edge_in_quad() = face;
	}
      nci += 4;
    }
  //nboundary += nci;
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::new_vertexs3d(HangContainer3d& hangset, 
				     const IntVector& vnew,
				     const IntSet& CellRefList)
{
  int nv1 = CellRefList.size();

  IntSetIt  cp=CellRefList.begin();
  for (int i=0; i<nv1; i++)
    {
      int f = co2n[*cp++];
      new_middle_vertex3d(vnew[i],f);
    }
  /* new vertexes on faces */  
  HangList<4>::const_iterator f = hangset.FaceCreating().begin();

  for (;f!=hangset.FaceCreating().end(); f++)
    {
      new_face_vertex3d(f->second.hanging(),f->first);
    }
  /* new vertexes on edges */
  HangList<2>::const_iterator p = hangset.Creating().begin();

  for (;p!=hangset.Creating().end(); p++)
    {
      new_edge_vertex3d(p->second.hanging(),p->first);
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::new_middle_vertex3d(int i, int f)
{
  int v0 = hexs[f].vertex(0);
  int v1 = hexs[f].vertex(1);
  int v2 = hexs[f].vertex(2);
  int v3 = hexs[f].vertex(3);
  int v4 = hexs[f].vertex(4);
  int v5 = hexs[f].vertex(5);
  int v6 = hexs[f].vertex(6);
  int v7 = hexs[f].vertex(7);

  Vertex3d w1, w2;

  w1.equ(0.25,vertexs3d[v0], 0.25,vertexs3d[v1],
	 0.25,vertexs3d[v2], 0.25,vertexs3d[v3]);
  w2.equ(0.25,vertexs3d[v4], 0.25,vertexs3d[v5],
	 0.25,vertexs3d[v6], 0.25,vertexs3d[v7]);

  vertexs3d[i].equ(0.5,w1,0.5,w2);
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::new_hexs(const HangContainer3d& hangset,
				  const IntVector& cnew, 
				  const IntVector& vnew, int nvold,
				  const IntSet& CellRefList)
{
  // neue zellen erzeugen
  // eintragen der "Vater-Vertexs" in den kindern

  int nci = 0;
  int ivm = 0;

  for (IntSetCIt  cp=CellRefList.begin();
       cp!=CellRefList.end(); cp++)
    {
      int father = co2n[*cp];

      assert(father!=-1);

      vector<int>&  qc = hexs[father].childs();
      qc.resize(8);

      int childlevel = hexs[father].level()+1;
      for(int ic=0; ic<8; ic++)
	{
	  int inold = cnew[nci+ic];
	  qc[ic] = inold;
	  hexs[inold].level()  = childlevel;
	  hexs[inold].father() = father;
	  hexs[inold].childs().resize(0);
	  hexs[inold].edges() = -1;
	}
      HexLaO.fill_corner_vertex_in_childs(hexs[father]);

      HexLaO.fill_middle_vertex_in_childs(hexs[father],vnew[ivm]);
      ivm++;
      nci += 8;

      // Edge Vertex -- linehanginfo schon ok (hanging) !
      int ive(-1);
      FaceVector faceglob;
      for (int i=0; i<6; i++)
	{
	  HexLaO.global_face_unsorted(faceglob,hex(father),i);
	  
	  ive = hangset.vertex_index(faceglob);

	  assert(ive>=0);

	  HexLaO.fill_face_vertex_in_childs(hexs[father],i,ive);
	}
      fixarray<2,int> lineglob;
      for (int i=0; i<12; i++)
	{
	  HexLaO.global_edge_unsorted(lineglob,hex(father),i);
	  
	  ive = hangset.vertex_index(lineglob);

	  if (ive<0)
	    {
	      cout << "Line " << lineglob << endl;
	      cout << "new vertex " << ive << endl;
	      cout << "father hex " << father << endl;
	    }

	  assert(ive>=0);

	  HexLaO.fill_edge_vertex_in_childs(hexs[father],i,ive);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::inner_vertex_newton3d(const IntVector& vnew,
					       const IntSet& CellRefList, 
					       const IntSet& adjustvertex)
{
  if (GetCurvedShapes().empty()) return;

  // baue lokalen set auf um spaeter nicht alle hex zu justieren
  IntSet Hexset;
  for (int i=0; i<Bquads.size(); i++)
    {
      int hi = Bquads[i].of_quad();
      Hexset.insert(hi);
    }

  fixarray<4,int> v;

  IntSetIt  cp=CellRefList.begin();

  for (int i=0; i<CellRefList.size(); i++)
    {
      int hi = co2n[*cp++];
      const Hex& h = hex(hi);

      if (Hexset.find(hi)==Hexset.end()) continue;


      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEU
      // edges
      for (int e=0;e<12;++e)
	{
	  int ev = HexLaO.edge_vertex(h,e);
	  if (adjustvertex.find(ev)!=adjustvertex.end()) continue;

	  fixarray<2,int> fe;
	  HexLaO.globalvertices_of_edge(h,fe,e);
	  vertexs3d[ev]  = vertexs3d[fe[0]];
	  vertexs3d[ev] += vertexs3d[fe[1]];
	  vertexs3d[ev]*=0.5;
	}
      // faces
      for (int f=0;f<6;++f)
	{
	  int fv = HexLaO.face_vertex(h,f);
	  if (adjustvertex.find(fv)!=adjustvertex.end()) continue;

	  fixarray<4,int> fe;
	  HexLaO.LoadEdgeVerticesOfFace(h,f,fe);
	  vertexs3d[fv].equ(0.25, vertexs3d[fe[0]],
			    0.25, vertexs3d[fe[1]],
			    0.25, vertexs3d[fe[2]],
			    0.25, vertexs3d[fe[3]]);
	}
      // middle
      int fv = HexLaO.middle_vertex(h);
      assert (adjustvertex.find(fv)==adjustvertex.end());
      fixarray<6,int> fe;
      HexLaO.LoadFaceVertices(h,fe);
      vertexs3d[fv]=0;
      for (int i=0;i<6;++i) vertexs3d[fv]+=vertexs3d[fe[i]];
      vertexs3d[fv]*=1./6.;
      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEU







      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ALT
//       for (int face=0; face<6; face++)
// 	{
// 	  HexLaO.LoadEdgeVerticesOfFace(h,face,v);
// 	  int mv = HexLaO.face_vertex(h,face);
// 	  if (adjustvertex.find(mv)!=adjustvertex.end()) continue;
// 	  new_face_vertex3d(mv,v);
// 	}
// //       fixarray<6,int> w;
// //       int mv = HexLaO.middle_vertex(h);
// //       HexLaO.LoadFaceVertices(h,w);
// //       new_vertex3d(mv,w);
      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ALT
}
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::basic_fill_neighbours3d()
{
  int n=0;
  FaceVector faceglob;
  for(int ic=0; ic<hexs.size(); ic++)
    {
      for(int i=0; i<6;i++)
	{
	  HexLaO.global_face_unsorted(faceglob,hex(ic),i);
	  HangList<4>::iterator p = QuadHang.find(faceglob);
	  
	  if(p!=QuadHang.end())
	    {
	      int cn = p->second.cneighbour();
	      int rn = p->second.rneighbour();

	      if( (cn==-1) && (!hexs[ic].sleep()) )
		{
		  n++;
		  p->second.cneighbour() = ic;
		}
	      else if(rn==-1)
		{
		  n++;
		  p->second.rneighbour() = ic;
		}
	    }
	 }
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::init_quad(BoundaryQuad& newquad)
{
  FaceVector v;
  IntSet newquad_set;
  newquad_set.insert(newquad[0]);
  newquad_set.insert(newquad[1]);
  newquad_set.insert(newquad[2]);
  newquad_set.insert(newquad[3]);
  for (int i=0; i<hexs.size(); i++)
    {
      for (int face=0; face<6; face++)
	{
	  HexLaO.global_face_unsorted(v,hex(i),face);
          IntSet v_set;
          v_set.insert(v[0]);
          v_set.insert(v[1]);
          v_set.insert(v[2]);
          v_set.insert(v[3]);
          if(newquad_set==v_set)
	    {
              HexLaO.global_face_unsorted(v,hex(i),face);
              newquad.vertex()=v;
              newquad.of_quad()      = i;
              newquad.edge_in_quad() = face;
              return;
	    }
	}
    }
  cerr << "BoundaryQuad not found !" << endl;
  cerr << newquad;
  abort();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::init_edges3d()
{
  if (withfaces)
    {
      FaceManager EM(edges,hexs,co2n,eo2n);
      EM.InitFaces();
      EM.SortHangings();
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::GetVertexesOfFace(fixarray<5,int>& v, int e) const
{
  const Edge& E = edge(e);
  const Hex* Q = &hex(E.master());

  int le = E.LocalMasterIndex();
  v[0] = (*Q)[le]; 
  v[1] = (*Q)[(le+1)%8]; 
  v[2] = (*Q)[(le+2)%8]; 
  v[3] = (*Q)[(le+3)%8]; 
  v[4] = -1;

  if(!Q->sleep()) 
    {
      if (E.slave()==-1) return;
      Q = &hex(E.slave());
      le = E.LocalSlaveIndex();
    }
  v[4] = HexLaO.face_vertex(*Q,le);
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::GetVertexesOfFace(fixarray<4,int>& v, int e) const
{
  const Edge& E = edge(e);
  const Hex*  Q = &hex(E.master());

  int le = E.LocalMasterIndex();
  v[0] = (*Q)[le]; 
  v[1] = (*Q)[(le+1)%8]; 
  v[2] = (*Q)[(le+2)%8]; 
  v[3] = (*Q)[(le+3)%8]; 
}

/*------------------------------------------------------*/

int    HierarchicalMesh3d::Vater(const int i) const
{
  return hex(i).father();
}


/*------------------------------------------------------*/

int HierarchicalMesh3d::nactivedescendants(int i)      const
{
  if (!hexs[i].sleep())
    return 1;
  int k=0;
  for (int j=0;j<hexs[i].nchilds();++j)
    k+=nactivedescendants(hexs[i].child(j));
  return k;
}

/*------------------------------------------------------*/

IntVector HierarchicalMesh3d::GetVertices(int c) const
{  
  IntVector v;
  for (int i=0;i<hexs[c].nvertexs();++i)
    v.push_back(hexs[c][i]);
  return v;
}

/*------------------------------------------------------*/

IntVector HierarchicalMesh3d::Nachkommen(const int i) const
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

IntVector HierarchicalMesh3d::Kinder     (const int i) const
{
  IntVector k = hex(i).childs();
  return k;
} 

/*------------------------------------------------------*/

IntVector HierarchicalMesh3d::Geschwister(const int i) const
{
  const Hex& q = hex(i);
  int father = q.father();
  if (father==-1)
    {
      IntVector n(1,i);
      return n;
    }
  return Kinder(father);
}


/*------------------------------------------------------*/


fixarray<4,int> HierarchicalMesh3d::ChildrenOfFace(int e) const
{
  int s  = edge(e).slave();
  int is = edge(e).LocalSlaveIndex();

  assert(s>=0);

  fixarray<4,int> f;
  for(int ii=0;ii<4;ii++)
    {
      int ic = hex(s).child(HexLaO.ChildsOfFace(is,ii));
      f[ii] = hex(ic).edge(HexLaO.ChildFace(is));
    }
  return f;
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::GetAwakeCells(set<int>& v) const
{
  for (int i=0; i<ncells(); i++)
    {
      if (!hexs[i].sleep()) v.insert(i);
    }
}


/*---------------------------------------------------*/

void HierarchicalMesh3d::GetAwakePatchs(set<int>& v) const
{
  for (int i=0; i<ncells(); i++)
    {
      if (!hexs[i].sleep()) 
	{
	  int f = hexs[i].father();
	  assert(f!=-1);
	  v.insert(f);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::ConstructQ2PatchMesh(IntVector& q2patchmesh) const
{
  typedef set<int>::iterator It;
  set<int> patche;
  GetAwakePatchs(patche);
  vector<set<int> > patch_on_level(nlevels());
  q2patchmesh.resize(0);
  for (It it = patche.begin();it!=patche.end();++it)
    patch_on_level[hexs[*it].level()].insert(*it);
  // grobgitterzellen kommen ins patchmesh
  for (It it=patch_on_level[0].begin();
       it!=patch_on_level[0].end();++it)
    q2patchmesh.push_back(*it);
  // der Rest wird eins groeber
  for (int l=1;l<nlevels();++l)
    {
      It it = patch_on_level[l].begin();
      while (it!=patch_on_level[l].end())
	{
	  int v = hexs[*it].father();
	  assert(v!=-1);
	  q2patchmesh.push_back(v);
	  IntVector nk = Nachkommen(v);
	  for (int i=0;i<nk.size();++i)
	    patch_on_level[hexs[nk[i]].level()].erase(nk[i]);
	  it = patch_on_level[l].begin();
	}
    }
}

/*---------------------------------------------------*/

IntVector HierarchicalMesh3d::ConstructQ4Patch(int c) const
{
  IntVector patch(125,-1);
  for(int i=0; i<125; i++)
  {
    // Vertex i steht an Position (x,y,z)
    int x = i%5;
    int y = (i%25)/5;
    int z = i/25;

    // Position von erstem Kind
    int fcx = x/3;
    int fcy = y/3;
    int fcz = z/3;
    // Index davon
    int fci = fcz*4+fcy*2+abs(fcx-fcy);

    // Position von Kind im Kind
    int scx = (x-2*fcx)/2;
    int scy = (y-2*fcy)/2;
    int scz = (z-2*fcz)/2;
    // Index davon
    int sci = scz*4+scy*2+abs(scx-scy);

    // Position des Vertex
    int vx = x-2*fcx-scx;
    int vy = y-2*fcy-scy;
    int vz = z-2*fcz-scz;
    // Index davon
    int vi = vz*4+vy*2+abs(vx-vy);

    patch[i] = hexs[hexs[hexs[c].child(fci)].child(sci)].vertex(vi);
  }
  return patch;
}

/*---------------------------------------------------*/

int HierarchicalMesh3d::NodeOnFace(int e) const
{
  // only for real hanging nodes
  int m  = edge(e).master();
  int im = edge(e).LocalMasterIndex();
  int s  = edge(e).slave();
  int is = edge(e).LocalSlaveIndex();

  if(s<0) 
    {
      assert(hex(im).sleep());
      return HexLaO.face_vertex(hex(m),im);
    }

  return HexLaO.face_vertex(hex(s),is);
}

/*------------------------------------------------------*/

int HierarchicalMesh3d::neighbour(int c, int le) const
{
  const Hex&  Q = hex(c);
  const Edge& E = edge(Q.edge(le));
  int m = E.master();
  int nq = m;
  if (m==c) nq = E.slave();
  return nq;
}

/*---------------------------------------------------*/

int HierarchicalMesh3d::neighbour_neighbour(int c, int le) const
{
  assert(le<6);
  int n = neighbour(c,le);
  assert(n>=0);

  int nn=0;
  for (nn=0;nn<6;++nn) if (c==neighbour(n,nn)) break;
  assert(nn<6);
  return nn;
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::FillAllBoundaryLines()
{
  cout << "HierarchicalMesh3d::FillAllBoundaryLines() not written" << endl;
  abort();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::WriteAll(const string& name) const
{
  ofstream out(name.c_str());

  out << dimension()  << " dimension" << endl;
  out << nnodes() << " vertexs"   << endl;
  out << mnlevels << " mnlevels"   << endl;

  cerr << "HierarchicalMesh3d::WriteAll()\n";
  cerr << "not written 3d\n";

  out.close();
  abort();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::write_inp(const string& name) const
{
  ofstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh3d::write_inp()\n";
      cerr << "cannot open file " << name << endl;
      abort();
    }
  
  int nt = ncells()+nbquads();
  file <<nnodes()<<" "<<nt<<" "<<0<<" "<<0<<" "<<0<<endl;

  
  for(int i=0;i<nnodes();i++) file <<i<<" "<<vertex3d(i) << " " << endl;

  for(int i=0;i<ncells();i++)
    {
      file << i<<" "<<0<<" hex "<<hex(i).vertex()<<endl;
    }
  for(int i=0;i<nbquads();i++)
    {
      file << i<<" "<< bquad(i).material() <<" quad "<<bquad(i).vertex()<<endl;
    }
}

/*---------------------------------------------------*/

pair<bool,tint> HierarchicalMesh3d::check_inp(const string& name)
{
  //  detect some errors in input-file... and more

  ifstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh3d::check_inp()\n";
      cerr << "cannot open file " << name << endl;
      abort();
    }

  bool first_one = 1;
  int  nv, nl, nq, nh, nt;
  int  n_unkonwn;
  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  Vertex3d c; int ind;
  for(int i=0;i<nv;i++)
    {
      file >> ind >> c;
    }

  nh = 0; nq = 0; nl = 0;
  fixarray<8,int> ih;
  fixarray<4,int> iq;
  fixarray<2,int> il;
  for(int i=0;i<nt;i++)
    {
      string name;
      string mat;
      int ii;
      file >> ii >> mat  >> name;
      if(name=="hex")
	{
	  file >> ih;
	  nh++;
	  if( (ih[0]==0)||(ih[1]==0)||(ih[2]==0)||(ih[3]==0) ||
	      (ih[4]==0)||(ih[5]==0)||(ih[6]==0)||(ih[7]==0))
	    {
	      first_one = 0;
	    }
	}
      else if(name=="quad")
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
  if(nt!=(nl+nq+nh))
    {
      cerr << "wrong number of cells: " << nt << endl;
      cerr << "lines quads hexs: " << nl << " " << nq << " " << nh << endl;
      abort();
    }

  file.close();

  return make_pair(first_one,make_triple(nl,nq,nh));
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::read_inp(const string& name)
{
  // check mesh.....

  pair<bool,tint> p = check_inp(name);
  bool first_one = p.first;
  tint n = p.second;

  int nl = n.first;
  int nq = n.second;
  int nh = n.third;

  ifstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh3d::read_inp()\n";
      cerr << "cannot open file " << name << endl;
      abort();
    }

  int  nv,  nt, n_unkonwn;

  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  cout << "3D Mesh:  " << nv << " nodes, ";
  cout << nl << " lines, ";
  cout << nq << " quads, ";
  cout << nh << " hexs" << endl;

  vertexs3d.reserve(nv);
  vertexs3d.resize(nv,Vertex3d());
  hexs.     reserve(nh);
  hexs.     resize(nh,Hex());
  Bquads.   reserve(nq);
  Bquads.   resize(nq);

  Vertex3d c;      int ind;
  for(int i=0;i<nv;i++)
    {
      file >> ind >> c;
      vertexs3d[i] = c;
    }
  fixarray<8,int> ihv;
  fixarray<4,int> iqv;
  fixarray<2,int> ilv;
  int ih = 0;
  int iq = 0;
  for(int i=0;i<nt;i++)
    {
      string name;	int unknown; string matstring;
      file >> unknown >> matstring  >> name;
      if(name=="hex")
	{
	  file >> ihv;
	  if(first_one) for(int iii=0;iii<8;iii++) ihv[iii]--;
	  hexs[ih++].vertex() = ihv;
	}
      else if(name=="quad")
	{
	  file >> iqv;
	  if(first_one) for(int iii=0;iii<4;iii++) iqv[iii]--;
	  
	  BoundaryQuad li;
	  li.material() = atoi(matstring.c_str());
	  li.vertex() = iqv;
	  init_quad(li);
	  Bquads[iq++]=li;
	}
    }
  init_edges3d();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::write(const string& bname) const
{
  write_gup(bname);
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::write_gup(const string& bname) const
{
  string name = bname;
  name += ".gup";

  ofstream out(name.c_str());

  out << dimension()  << " dimension" << endl;
  out << nnodes() << " vertexs"   << endl;

  for (int i=0; i<nnodes(); i++)
    {
      out << " " << vertex3d(i) << endl;
    }
  out << hexs.size() << " hexs" << endl;
  for (int i=0; i<hexs.size(); i++)
    {
      out << hex(i);
    }
  out << QuadHang << endl;
  out << LineHang << endl;
  out << Bquads.size() << " boundaryquads" << endl;
  for (int i=0; i<Bquads.size(); i++)
    {
      out << Bquads[i].material() << " " << Bquads[i] << endl;
    }
  out << endl << edges.size() << " edges" << endl;
  for (int i=0; i<edges.size(); i++)
    {
      out << " " << edges[i];
    }
  out.close();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::write_gip(const string& bname) const
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
      vertex3d(i).BinWrite(out);
    }

  int nhexs = hexs.size();
  out.write(reinterpret_cast<const char*>(&nhexs),sizeInt);
 
  for (int i=0; i<hexs.size(); i++)
    {
      hex(i).BinWrite(out);
    }

  QuadHang.BinWrite(out);
  LineHang.BinWrite(out);
  int nbquads=Bquads.size();
  out.write(reinterpret_cast<const char*>(&nbquads),sizeInt);
  for (int i=0; i<Bquads.size(); i++)
    {
      int mat = Bquads[i].material();
      out.write(reinterpret_cast<const char*>(&mat),sizeInt);
      Bquads[i].BinWrite(out);
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

void HierarchicalMesh3d::read_gup(const string& name)
{
  string symbol;

  ifstream file(name.c_str());

  assert(file.is_open());

  Bquads.clear();
  edges.clear();
  QuadHang.clear();
  LineHang.clear();

  int n;
  int dim;
  file >> dim >> symbol >> n >> symbol;

  assert(dim==3);

  if( _i_showoutput ){
    cout << "Mesh 3d  : ";
  }
  vertexs3d.reserve(n);
  vertexs3d.resize(n);

  assert(symbol=="vertexs");

  for (int i=0; i<n; i++)
    {
      file >> vertexs3d[i];
    }
  if( _i_showoutput ){
    cout << n << " nodes, ";
  }

  file >> n >> symbol;

  assert(symbol=="hexs");

  hexs.reserve(n);
  hexs.resize(n);
  for (int i=0; i<hexs.size(); i++)
    {
      file >> hexs[i];
    }
  if( _i_showoutput ){
    cout <<  n << " hexs, ";
  }
  file >> QuadHang;
  if( _i_showoutput ){
    cout << QuadHang.size() << " quadhangs, ";
  }
  file >> LineHang;
  if( _i_showoutput ){
    cout << LineHang.size() << " linehangs, ";
  }
  
  file >> n >> symbol;
  int number = 0;
  assert(symbol=="boundaryquads");

  BoundaryQuad bol;
  for (int i=0; i<n; i++)
    {
      file >> bol.material() >> bol;
      Bquads.push_back(bol);
    }
  number = n;

  if( _i_showoutput ){
    cout << number << " boundaryquads, ";
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
  if (edges.size()==0) init_edges3d();
  post_refine3d();
  file.close();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::read_gip (const string& bname)
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
      cerr << "HierarchicalMesh3d::read_gip(): error in file "<< name <<endl;
      abort();
    }

  Bquads.clear();
  edges.clear();
  QuadHang.clear();
  LineHang.clear();

  int n, dim, sizeInt=sizeof(int);
  file.read(reinterpret_cast<char*>(&dim),sizeInt);
  file.read(reinterpret_cast<char*>(&n),sizeInt);

  assert(dim==3);

  if( _i_showoutput ){
    cout << "Mesh 3d  : ";
  }
  vertexs3d.reserve(n);
  vertexs3d.resize(n);

  for (int i=0; i<n; i++)
    {
      vertexs3d[i].BinRead(file);
    }
  if( _i_showoutput ){
    cout << n << " nodes, ";
  }
  file.read(reinterpret_cast<char*>(&n),sizeInt);
  hexs.reserve(n);
  hexs.resize(n);
  for (int i=0; i<hexs.size(); i++)
    {
      hexs[i].BinRead(file);
    }
  if( _i_showoutput ){
    cout <<  n << " hexs, ";
  }
  QuadHang.BinRead(file);
  if( _i_showoutput ){
    cout << QuadHang.size() << " quadhangs, ";
  }
  LineHang.BinRead(file);
  if( _i_showoutput ){
    cout << LineHang.size() << " linehangs, ";
  }
  file.read(reinterpret_cast<char*>(&n),sizeInt);
  int number = 0;
  BoundaryQuad bol;
  for (int i=0; i<n; i++)
    {
      file.read(reinterpret_cast<char*>(&bol.material()),sizeInt);
      bol.BinRead(file);
      Bquads.push_back(bol);
    }
  number = n;
  if( _i_showoutput ){
    cout << number << " boundaryquads, ";
  }
  file.read(reinterpret_cast<char*>(&n),sizeInt);
  if( _i_showoutput ){
    cout << n << " edges" << endl;
  }
  Edge e;
  for (int i=0; i<n; i++)
    {
      e.BinRead(file);
      edges.push_back(e);
    }
  if (edges.size()==0) init_edges3d();
  post_refine3d();
  file.close();
}

/*---------------------------------------------------*/

int HierarchicalMesh3d::regular_grid3d_one(IntSet& celllist,IntVector& coarsesub,
					 const IntSet& CellRefList,
					 const IntSet& CellCoarseList)
{
  /* detects jump over two levels across LineHangs */

  int n = 0;
  HangList<4>::const_iterator  hp;
 
  for (hp = QuadHang.begin(); hp!=QuadHang.end(); ++hp)
  {
    int cr = hp->second.rneighbour();
    int cn = hp->second.cneighbour();

    assert(cr>=0);

    if(cn!=-1)
      {
	if(hex(cr).childs().size()==0) continue;

	FaceVector  f;
	HexLaO.childs_of_global_face(f,hex(cr),hp->first);

	for(unsigned i=0;i<f.size();++i)
	  {
	    int c = f[i];

	    if ( (CellRefList.find(c)!=CellRefList.end()) ||
		 (hexs[c].sleep()) )
	      {
		if (CellCoarseList.find(cn)!=CellCoarseList.end())
		  {
		    coarsesub.push_back(cn);
		    break;
		  }
		else
		  {
		    assert(!hex(cn).sleep());

		    pair<IntSetIt,bool> p = celllist.insert(cn);
		    if (p.second)  
		      { 
			n++;
			break; 
		      }
		    //celllist.push_back(cn);
		    //break; 
		  }
	      }
	  }
      }
  }
  return n+coarsesub.size();
  //unique(coarsesub.begin(),coarsesub.end());
  //unique(celllist.begin(),celllist.end());
  //return celllist.size()+coarsesub.size();
}

/*---------------------------------------------------*/

int HierarchicalMesh3d::regular_grid3d_one(IntVector& celllist,IntVector& coarsesub,
					 const IntSet& CellRefList,
					 const IntSet& CellCoarseList)
{
  IntSet h;
  Vec2Set(h,celllist);
  int r = regular_grid3d_one(h,coarsesub,CellRefList,CellCoarseList);
  Set2Vec(celllist,h);
  return r;

  /* detects jump over two levels across LineHangs */

  //int n = 0;
  HangList<4>::const_iterator  hp;
 
  for (hp = QuadHang.begin(); hp!=QuadHang.end(); ++hp)
  {
    int cr = hp->second.rneighbour();
    int cn = hp->second.cneighbour();
    assert(cr>=0);

    if(cn!=-1)
      {
	if(hex(cr).childs().size()==0) continue;

	FaceVector  f;
	HexLaO.childs_of_global_face(f,hex(cr),hp->first);

	for(unsigned i=0;i<f.size();++i)
	  {
	    int c = f[i];

	    if ( (CellRefList.find(c)!=CellRefList.end()) ||
		 (hexs[c].sleep()) )
	      {
		if (CellCoarseList.find(cn)!=CellCoarseList.end())
		  {
		    coarsesub.push_back(cn);
		    break;
		  }
		else
		  {
		    assert(!hex(cn).sleep());
		    int number =0;
		    for (int kk=0; kk<celllist.size(); kk++)
		      {
			if (celllist[kk]==cn) number++;
		      }
		    if (!number)
		      {
			celllist.push_back(cn);
			//n++;
			break; 
		      }
		  }
	      }
	  }
      }
  }
  return celllist.size()+coarsesub.size();
}

/*---------------------------------------------------*/

int HierarchicalMesh3d::regular_grid3d_two(IntVector& celllist,
					 const IntSet& CellRefList)
{
  /* detects more than 4 HangFaces on one Hex */

  IntVector nh(hexs.size());
  for (HangList<4>::const_iterator p = QuadHang.begin();
       p!=QuadHang.end(); p++)
    {
      int i = p->second.cneighbour();
      if (i<0)                                    continue;
      if (hex(i).sleep())                         continue;
      if (CellRefList.find(i)!=CellRefList.end()) continue;

      nh[i]++;
    }
  for (int i=0; i<hexs.size(); i++)
    {
      if (nh[i]>4)
	{
	  //pair<IntSetIt,bool> pp = celllist.insert(i);
	  //if(pp.second)  nto++;
	  celllist.push_back(i);
	  //nto++;
	}
    }
  return celllist.size();
  //return nto;
}


/*---------------------------------------------------*/

void HierarchicalMesh3d::GetMinMaxLevels(IntVector& maxi, IntVector& mini,
					 const IntSet& CellRef) const
{
  // set maximal levels for vertices
  //
  maxi.resize(nnodes());
  mini.resize(nnodes()); mini = 1000;
  for (int i=0; i<hexs.size(); i++)
    {
      const Hex& q = hex(i);
      if(q.sleep()) continue;

      int lev = q.level();
      if (CellRef.find(i)!=CellRef.end()) lev++;
      for (int j=0; j<8; j++)
	{
	  int k = q[j];
	  maxi[k] = Gascoigne::max_int( maxi[k], lev);
	  mini[k] = Gascoigne::min_int( mini[k], lev);
	}
    }
}

/*---------------------------------------------------*/

int HierarchicalMesh3d::regular_grid3d_three_refine(IntSet& CellRef) const
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
  for (int i=0; i<hexs.size(); i++)
    {
      const Hex& q = hex(i);
      if (q.sleep()) continue;
      int lev = q.level();
      for (int j=0; j<8; j++)
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
  return ref;
}

/*---------------------------------------------------*/

int HierarchicalMesh3d::regular_grid3d_three_coarse(IntSet& CellRef,
						    IntSet& CellCoarse) const
{
  int maxl = 0;
  vector<IntSet> LevelCellCoarse;
  {
    IntSet::const_iterator p = CellCoarse.begin();
    for ( ; p!= CellCoarse.end(); p++)
      {
	const Hex& q = hex(*p);
	int lev = q.level();
	maxl = Gascoigne::max_int(maxl,lev);
      }
    LevelCellCoarse.resize(maxl+1);
    p = CellCoarse.begin();
    for ( ; p!= CellCoarse.end(); p++)
      {
	const Hex& q = hex(*p);
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
	  const Hex& q = hex(*p);
	  int lev = q.level();
	  for (int j=0; j<8; j++)
	    {
	      int k = q[j];
	      if (lev+1<maxlevel[k])
		{
		  coarsesub.insert(*p);
		  continue;
		}
	    }
	}
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
  return coarse;
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::global_coarse3d()
{
  IntSet RefList,CoarseList;

  for (int i=0; i<hexs.size(); i++)
    {
      const Hex& q = hexs[i];
      if ( q.sleep()) continue;
      if (!q.level()) continue;
      CoarseList.insert(q.father());
    }

  HangContainer3d hangset(LineHang,QuadHang);

  UpdateHangs(hangset,RefList,CoarseList);

  basic_refine3d(hangset,RefList,CoarseList);

  post_refine3d();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::LoadFathers3d(IntVector& v) const
{
  IntSet fathers;
  
  for (int i=0; i<v.size(); i++)
    {
      if (hex(v[i]).level()==0) continue;
      int f = hex(v[i]).father();

      assert(f>=0) ;

      fathers.insert(f);
    }
  Set2Vec(v,fathers);
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::FillVertexLevels(IntVector& dst) const
{
  dst.resize(nnodes(),100);
  for (int i=0; i<ncells(); i++)
    {
      const Hex& Q = hexs[i];
      if (Q.sleep()) continue;

      int level = Q.level();
      for (int j=0; j<8; j++)
	{
	  int k = Q[j];
	  dst[k] = Gascoigne::min_int(dst[k],level);
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::RefineCoarseNodes
(IntSet& dst, const IntVector& refnodes, const IntVector& vertexlevel) const
{
  IntSet h;

  Vec2Set(h,refnodes);
  
  dst.clear();
  for (int i=0; i<ncells(); i++)
    {
      const Hex& Q = hexs[i];
      if (Q.sleep()) continue;

      int f = Q.father();
      if (f<0)  continue;

      const Hex& QF = hexs[f];
	   
      for (int j=0; j<8; j++)
	{
	  int k = Q[j];
	  if (h.find(k)==h.end()) continue;
	  
	  int minlevel = vertexlevel[QF[0]];
	  for (int v=1; v<8; v++)
	    {
	      minlevel = Gascoigne::min_int(minlevel,vertexlevel[QF[v]]);
	    }
	  for (int v=0; v<8; v++)
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

void HierarchicalMesh3d::VertexToCells(IntVector& dst, const IntSet& src, 
				       const IntVector& vertexlevel) const
{
  for (int i=0; i<ncells(); i++)
    {
      const Hex& Q = hexs[i];
      if (Q.sleep()) continue;
      int level = Q.level();
      for (int j=0; j<8; j++)
	{
	  int k = Q[j];
	  if (vertexlevel[k]==level)
	    {
	      if (src.find(k)!=src.end())
		{
		  dst.push_back(i);
		}
	    }
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::VertexToCellsCoarsening(IntVector& dst, const IntSet& src, 
						 const IntVector& vertexlevel) const
{
  for (int i=0; i<ncells(); i++)
    {
      const Hex& Q = hexs[i];
      if (Q.sleep()) continue;
      int level = Q.level();
      int count = 0;
      for (int j=0; j<8; j++)
	{
	  int k = Q[j];
	  if (vertexlevel[k]==level)
	    {
	      if (src.find(k)!=src.end())
		{
		  count++;
		}
	      else break;
	    }
	}
      if (count==8) dst.push_back(i);
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::recursive_childs(int q, IntVector& ref, int d) const
{
  if (d>0)
    {
      const Hex& Q = hex(q);
      assert(Q.sleep());
      for (int j=0; j<8; j++)
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

void HierarchicalMesh3d::patch_refine(IntVector& cell_ref, 
				      IntVector& cell_coarse)
{
  for (int i=0; i<pdepth; i++)
    {
      LoadFathers3d(cell_ref);
      LoadFathers3d(cell_coarse);
    }

  CoarseHierarchicalMesh3d CM(*this);

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

void HierarchicalMesh3d::refine
(const IntVector& cell_ref, const IntVector& cell_coarse)
{
  IntSet CellRefList, CellCoarseList;
  _refine3d(CellRefList,CellCoarseList,cell_ref,cell_coarse);
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::_refine3d
(IntSet& CellRefList, IntSet& CellCoarseList, 
 const IntVector& cell_ref, const IntVector& cell_coarse)
{
  prepare3d(cell_ref,cell_coarse, CellRefList, CellCoarseList);

  while (regular_grid3d_three_refine(CellRefList)) {}
  while (regular_grid3d_three_coarse(CellRefList,CellCoarseList)) {}

  HangContainer3d hangset(LineHang,QuadHang);

  UpdateHangs(hangset,CellRefList,CellCoarseList);

  basic_refine3d(hangset,CellRefList,CellCoarseList);

  post_refine3d();
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::delete_vertexs3d(const IntVector& vo2n)
{
  for(unsigned oi=0; oi<vo2n.size(); ++oi)
    {
      int ni = vo2n[oi];
      if(ni>=0)  
	{
	  vertexs3d[ni] = vertexs3d[oi];
	}
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::post_refine3d()
{
  //check_mesh3d();
  mnlevels = 0;
  for(int i=0;i<hexs.size();i++)
    {
      mnlevels = Gascoigne::max_int(mnlevels,hexs[i].level());
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::new_edge_vertex3d(int nv, const EdgeVector& v)
{
  vertexs3d[nv].equ(0.5, vertexs3d[v[0]], 0.5, vertexs3d[v[1]]);
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::new_face_vertex3d(int nv, const FaceVector& v)
{
  vertexs3d[nv].equ(0.25, vertexs3d[v[0]], 0.25, vertexs3d[v[1]],
		  0.25, vertexs3d[v[2]], 0.25, vertexs3d[v[3]]);
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::new_vertex3d(int nv, const fixarray<6,int>& v)
{
  cout << "@ " << vertexs3d[nv] << "\t";
  double d = 1./6.;
  vertexs3d[nv].equ(d, vertexs3d[v[0]], d, vertexs3d[v[1]],
		    d, vertexs3d[v[2]], d, vertexs3d[v[3]],
		    d, vertexs3d[v[4]], d, vertexs3d[v[5]]);
  cout << " -> " << vertexs3d[nv] << "\n";
}

/*---------------------------------------------------*/

void HierarchicalMesh3d::check_mesh3d() const
{

  // check quads

  int cmin = 0, cmax = hexs.size(); 
  int vmin = 0, vmax = nnodes(); 

  for(vector<Hex>::const_iterator p = hexs.begin();
      p!=hexs.end();p++)
    {
      const Hex& q = *p;

      // check vertex id
      for(int i=0; i< q.nvertexs(); i++)
	{
	  if((q.vertex(i)<vmin)||(q.vertex(i)>vmax))
	    {
	      cerr << "Vertex invalid in Cell: " << *p << " : ";
	      cerr << q.vertex(i) << endl;
	      exit(1);
	    }
	}
      // check child id
      for(int i=0; i< q.nchilds(); i++)
	{
	  if((q.child(i)<cmin)||(q.child(i)>cmax))
	    {
	      cerr << "Chid invalid in Cell: " << *p << " : ";
	      cerr << q.child(i) << endl; 
	      exit(1);
	    }
	}
    }

  // check quadhang
  for (HangList<4>::const_iterator  hp = QuadHang.begin(); hp!=QuadHang.end(); ++hp)
  {
    int cr = hp->second.rneighbour();
    if(cr==-1)
      {
	cerr << "Refine Neighbour invalid in hang: ";
	//cerr << *hp << endl;
	exit(1);
      }
  }
}
}

/*---------------------------------------------------*/

