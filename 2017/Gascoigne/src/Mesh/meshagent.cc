/**
*
* Copyright (C) 2004, 2005, 2006, 2008, 2009 by the Gascoigne 3D authors
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


#include  "meshagent.h"
#include  "visualization.h"
#include  "filescanner.h"
#include  "gascoignemeshconstructor.h"
#include  "stringutil.h"
#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HASHMAP std::tr1::unordered_map
#else
#include  <ext/hash_map>
#define HASHMAP __gnu_cxx::hash_map
#endif

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  MeshAgent::MeshAgent() : MeshAgentInterface(), _goc2nc(false), HMP(NULL), GMG(NULL)
{
}

/*-----------------------------------------*/

MeshAgent::~MeshAgent()
{
  if (HMP!=NULL) { delete HMP; HMP=NULL;}
  if (GMG!=NULL) { delete GMG; GMG=NULL;}
}

/*-----------------------------------------*/

void MeshAgent::ReInit()
{
  GMG->ReInit(GetDimension(),HMP->nlevels()-HMP->patchdepth());

  GascoigneMeshConstructor MGM(HMP,GMG);
  MGM.BasicInit();
  _celll2g = MGM.Celll2g();
  _cellg2l = MGM.Cellg2l();
  if(HMP->patchdepth()>=2)
  {
    BuildQ4PatchList(MGM.Patchl2g());
    assert(_q4patch.size()==_q4toq2.size());
    PatchIndexHandler &PIH = GMesh(0)->GetPatchIndexHandler();
    nvector<IntVector> &q4patch2cell = PIH.GetAllQ4Patch2Cell();
    q4patch2cell.clear();
    for(int i=0; i<_q4patch.size(); i++)
    {
      IntVector q2p = _q4toq2[i];
      IntVector q4p2c(0);
      for(int j=0; j<q2p.size(); j++)
      {
        int p = q2p[j];
        IntVector cells = PIH.GetPatch2Cell(p);
        q4p2c.insert(q4p2c.end(),cells.begin(),cells.end());
      }
      q4patch2cell.push_back(q4p2c);
    }
    PIH.GetIndexQ4() = _q4patch;
  }

  if(_goc2nc)
  {
      _co2n.clear();
      for(int i =0 ; i < _cl2g.size(); i++)
      {
	  int cn_i = _cl2g[i];//HM Num nach alter nummerierung
	  cn_i = HMP->Cello2n(cn_i);//HM Num nach neuer nummerierung
	  if(cn_i <0)
	  {
	      //Hier wurde vergroebert
	  }
	  else
	  {
	      set<int> kinder;
	      // Zuordnung alte GM Nummern zu neuen
	      if(HMP->sleep(cn_i))
	      {
		  //Zelle verfeinert
		  for(int j = 0; j <HMP->nchilds(cn_i); j++)
		  {
		      kinder.insert(_cellg2l[HMP->child(cn_i,j)]);
		  }
	      }
	      else
	      {
		  kinder.insert(_cellg2l[cn_i]);
	      }
	      _co2n[i] = kinder;
	  }
      }
      

      // Die Var fuer die alten Werteumschreiben
      _cl2g = _celll2g;
      _cg2l = _cellg2l;
      
      //den fathers Vektor neu fuellen
      _fathers.clear();
      _fathers.resize(_cl2g.size());
      for(int i = 0; i < _fathers.size(); i++)
      {
	  _fathers[i] = HMP->Vater(_cl2g[i]);
      }
  }

  AssemblePeriodicBoundaries();
}

/*-----------------------------------------*/

void MeshAgent::BuildQ4PatchList(const IntVector &patchl2g)
{
  _q4patch.resize(0);
  _q4toq2.resize(0);
  IntSet q2patch,q4patchset;
  HMP->GetAwakePatchs(q2patch);
  for(IntSet::const_iterator p=q2patch.begin(); p!=q2patch.end(); p++)
  {
    assert(*p!=-1);
    int vater = HMP->Vater(*p);
    assert(vater!=-1);
    q4patchset.insert(vater);
  }
  for(IntSet::const_iterator p=q4patchset.begin(); p!=q4patchset.end(); p++)
  {
    _q4patch.push_back(HMP->ConstructQ4Patch(*p));
  }
  if(patchl2g.size()!=q2patch.size())
  {
    cerr << "MeshAgent::BuildQ4PatchList: patchl2g must be same size as q2patch!!!" << endl;
    abort();
  }
  HASHMAP<int,int> patchg2l;
  for(int i=0; i<patchl2g.size(); i++)
  {
    //Edit patchg2l[patchl2g[i]] = i;
    patchg2l.insert(make_pair<int,int>(patchl2g[i],i));
  }
  IntVector perm(4);
  perm[0]=0;
  perm[1]=1;
  perm[2]=3;
  perm[3]=2;
  if(HMP->dimension()==3)
  {
    perm.resize(8);
    perm[4]=4;
    perm[5]=5;
    perm[6]=7;
    perm[7]=6;
  }
  _q4toq2.resize(q4patchset.size());
  int q4_l = 0;
  for(IntSet::const_iterator p=q4patchset.begin(); p!=q4patchset.end(); p++,q4_l++)
  {
    const IntVector &K = HMP->Kinder(*p);
    assert(K.size()==perm.size());
    for(int i=0; i<K.size(); i++)
    {
      assert(patchg2l.find(K[perm[i]])!=patchg2l.end());
      _q4toq2[q4_l].push_back(patchg2l[K[perm[i]]]);
    }
  }
}

/*-----------------------------------------*/

void MeshAgent::AssemblePeriodicBoundaries()
{
  /*-------------------------------------------------------
  | Erstellt eine Zuordnung der Knoten der periodischen
  | Raender zueinander.
  | Wird fuer die Modifikation der Systemmatrix und der
  | rechten Seite benoetigt; legt fest, welche Zeilen und
  | Spalten bzw. Eintraege miteinander kombiniert werden
  | muessen.
  -------------------------------------------------------*/
  int i_dim = GetDimension();

  if (_periodicCols.size() != 0)
  {
    int n = GMG->nlevels();

    for (int i = 0; i < n; i++)
    {
      MeshInterface* p_mesh = GMG->GetGascoigneMesh(i);
      GascoigneMesh* GMP    = dynamic_cast<GascoigneMesh*>(p_mesh);
      assert(GMP);
 
      //TODO: das soll eigentlich zu Solver gehoeren, damit die anderen Member darauf zugreifen koennen
      map<int,map<int,int> > mm_PeriodicPairs;
 
      assert(_periodicCols.size()%2 == 0);
 
      for(IntVector::const_iterator p_col   = _periodicCols.begin(); p_col!=_periodicCols.end(); p_col++)
      {
        map<int,int> m_PeriodicPairsCol1, m_PeriodicPairsCol2;
        int col_v = *p_col;
        int col_w = *(++p_col);
        
        if (_periodicMaps.find(col_v) == _periodicMaps.end())
        {
          cerr << "Periodic mapping not found: " << col_v << " --> " << col_w << endl;
          abort();
        }
        if (_periodicMaps[col_v].find(col_w) == _periodicMaps[col_v].end())
        {
          cerr << "Periodic mapping not found: " << col_v << " --> " << col_w << endl;
          abort();
        }
  
        IntVector iv_FirstBoundary  =  *(GMP->VertexOnBoundary(col_v));
        IntVector iv_SecondBoundary =  *(GMP->VertexOnBoundary(col_w));
  
        double max_diff_bestfit = 0.;
   
        // process each node on the first boundary
        for (IntVector::const_iterator p_first = iv_FirstBoundary.begin(); p_first != iv_FirstBoundary.end();p_first++)
        {
          if (i_dim == 2)
          {
            Vertex2d first = GMP->vertex2d(*p_first);
            Vertex2d otherside;
            _periodicMaps[col_v][col_w]->transformCoords(otherside, first);
            int bestfit;
            double diff, diff_bestfit;
            Vertex2d second;
            IntVector::const_iterator p_second = iv_SecondBoundary.begin();
            second = GMP->vertex2d(*p_second);
            bestfit = *p_second;
            diff_bestfit = (second.x()-otherside.x())*(second.x()-otherside.x())
              + (second.y()-otherside.y())*(second.y()-otherside.y());
   
            // find the best fit on the other boundary
            for (; p_second != iv_SecondBoundary.end(); p_second++)
            {
              second = GMP->vertex2d(*p_second);
              diff = (second.x()-otherside.x())*(second.x()-otherside.x())
                   + (second.y()-otherside.y())*(second.y()-otherside.y());
              if (diff < diff_bestfit)
              {
                bestfit = *p_second;
                diff_bestfit = diff;
              }
            }
            if (diff_bestfit > max_diff_bestfit) { max_diff_bestfit = diff_bestfit;}
   
            // create a map for fast access
            m_PeriodicPairsCol1[*p_first] = bestfit;
            m_PeriodicPairsCol2[bestfit]  = *p_first;
          }
          else if (i_dim == 3)
          {
            Vertex3d first = GMP->vertex3d(*p_first);
            Vertex3d otherside;
            _periodicMaps[col_v][col_w]->transformCoords(otherside, first);
            int bestfit;
            double diff, diff_bestfit;
            Vertex3d second;
            IntVector::const_iterator p_second = iv_SecondBoundary.begin();
            second = GMP->vertex3d(*p_second);
            bestfit = *p_second;
            diff_bestfit = (second.x()-otherside.x())*(second.x()-otherside.x())
                         + (second.y()-otherside.y())*(second.y()-otherside.y())
                         + (second.z()-otherside.z())*(second.z()-otherside.z());
 
            // find the best fit on the other boundary
            for (; p_second != iv_SecondBoundary.end(); p_second++)
            {
              second = GMP->vertex3d(*p_second);
              diff = (second.x()-otherside.x())*(second.x()-otherside.x())
                   + (second.y()-otherside.y())*(second.y()-otherside.y())
                   + (second.z()-otherside.z())*(second.z()-otherside.z());
              if (diff < diff_bestfit)
              {
                bestfit = *p_second;
                diff_bestfit = diff;
              }
            }
            if (diff_bestfit > max_diff_bestfit) { max_diff_bestfit = diff_bestfit;}
   
            // create a map for fast access
            m_PeriodicPairsCol1[*p_first] = bestfit;
            m_PeriodicPairsCol2[bestfit]  = *p_first;
          }
          else
          {
            std::cerr << "AssemblePeriodicBoundaries funktioniert nur fuer Dimension 2 und 3.\n";
            abort();
          }
         
          if (max_diff_bestfit > 1e-10 && _periodicMaps[col_v][col_w]->GetName() != "StdPeriodicMapping")
          {
            std::cerr << "Distance between boundaries " << col_v << " and " << col_w << " too large. Check your PeriodicMapping." << std::endl;
            abort();
          }
          
          mm_PeriodicPairs[col_v] = m_PeriodicPairsCol1;
          mm_PeriodicPairs[col_w] = m_PeriodicPairsCol2;
    
        }
      }

      GMP->GetBoundaryIndexHandler().SetPeriodicPairs(mm_PeriodicPairs);
    }
  }
}

/*-----------------------------------------*/

void MeshAgent::BasicInit(const ParamFile* paramfile)
{
  assert(HMP==NULL);
  int dim = 0;

  {
    DataFormatHandler DFH;
    DFH.insert("dimension",&dim);
    //um die zuordnung alte GMNr. -> GMNr. an/abzuschalten
    DFH.insert("cellnumtrans",&_goc2nc,false);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(paramfile,"Mesh");
  }
  {
    DataFormatHandler DFH;
    DFH.insert("periodic",&_periodicCols);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(paramfile,"BoundaryManager");
  }

  if (dim==2)
    {
      HMP = new HierarchicalMesh2d;
      for(map<int,BoundaryFunction<2>* >::const_iterator p=_curved2d.begin();p!=_curved2d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else if (dim==3)
    {
      HMP = new HierarchicalMesh3d;
      for(map<int,BoundaryFunction<3>* >::const_iterator p=_curved3d.begin();p!=_curved3d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else
    {
      cout << "dimension of Mesh ? " << dim << endl;
    }
  assert(HMP);
  HMP->BasicInit(paramfile);
  
  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::BasicInit(const string& gridname, int dim, int patchdepth, int epatcher, bool goc2nc)
{
  assert(HMP==NULL);
  _goc2nc = goc2nc;
  if (dim==2)
    {
      HMP = new HierarchicalMesh2d;
      for(map<int,BoundaryFunction<2>* >::const_iterator p=_curved2d.begin();p!=_curved2d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else if (dim==3)
    {
      HMP = new HierarchicalMesh3d;
      for(map<int,BoundaryFunction<3>* >::const_iterator p=_curved3d.begin();p!=_curved3d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else
    {
      cout << "dimension of Mesh ? " << dim << endl;
    }
  assert(HMP);
  HMP->SetParameters(gridname,patchdepth,epatcher);
  
  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::read_gup(const string& fname)
{
  assert(HMP);
  HMP->read_gup(fname);
  ClearCl2g();
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::read_gip(const string& fname)
{
  assert(HMP);
  HMP->read_gip(fname);
  ClearCl2g();
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::write_gup(const string& fname) const
{
  assert(HMP);
  HMP->write_gup(fname);
}

/*-----------------------------------------*/

void MeshAgent::write_gip(const string& fname) const
{
  assert(HMP);
  HMP->write_gip(fname);
}

/*-----------------------------------------*/

void MeshAgent::write_inp(const string& fname) const
{
  assert(HMP);
  HMP->write_inp(fname);
}

/*-----------------------------------------*/

void MeshAgent::global_patch_coarsen(int n)
{
  assert(HMP);
  HMP->global_patch_coarsen(n);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::global_refine(int n)
{
  assert(HMP);
  HMP->global_refine(n);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::random_patch_coarsen(double p, int n)
{
  assert(HMP);
  HMP->random_patch_coarsen(p,n);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::random_patch_refine(double p, int n)
{
  assert(HMP);
  HMP->random_patch_refine(p,n);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::refine_nodes(IntVector& refnodes)
{
  IntVector coarsenodes(0);
  refine_nodes(refnodes,coarsenodes);
}

/*-----------------------------------------*/

void MeshAgent::refine_nodes(IntVector& refnodes, IntVector& coarsenodes)
{
  assert(HMP);
  HMP->vertex_patch_refine(refnodes,coarsenodes);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::refine_cells(IntVector& ref)
{
  IntVector refnodes;
  
  for (int i=0; i<ref.size(); i++)
    {
      int cell = ref[i];
      for (int j=0; j<HMP->nodes_per_cell(cell); j++)
        {
          refnodes.push_back(HMP->vertex_of_cell(cell,j));
        }
    }
  refine_nodes(refnodes);
}

/*----------------------------------------*/

inline const set<int> MeshAgent::Cello2n(int i)const
{
    map<int,set<int> >::const_iterator p = _co2n.find(i);
    if(p == _co2n.end())
    {
	return set<int>();
    }
    else
    {
	return p->second;
    }
}

/*----------------------------------------*/

inline const int MeshAgent::Cello2nFather(int i)const
{
    assert(_co2n.find(i)==_co2n.end());
    //Umrechnung alte HM nummer in neue GM nummer
    return _cg2l.find(HMP->Cello2n(_fathers[i]))->second;
}


}

#undef HASHMAP
