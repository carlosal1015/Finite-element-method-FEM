/**
*
* Copyright (C) 2005, 2006, 2007 by the Gascoigne 3D authors
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


#include "backup.h"
#include "meshinterpolator.h"
#include "q12d.h"
#include "q13d.h"
#include "q22d.h"
#include "q23d.h"
#include "domainrighthandside.h"
#include <iterator>

using namespace std;

namespace Gascoigne
{
/**********************************************************/

class ProjectionRightHandSide : public DomainRightHandSide
{
  protected:

    int _ncomp;
    mutable FemFunction __U;

  public:

    ProjectionRightHandSide(int ncomp) : DomainRightHandSide(), _ncomp(ncomp) { }
    ~ProjectionRightHandSide() { }

    int GetNcomp() const { return _ncomp; }
    std::string GetName() const { return "ProjectionRightHandSide"; }

    void SetFemData(FemData& q) const
    {
      assert(q.count("U")==1);
      __U = q["U"];
    }

    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
    {
      for (int i=0; i<_ncomp; i++)
        b[i] += __U[i].m() * N.m();
    }
    void operator()(VectorIterator b, const TestFunction& N, const Vertex3d& v) const 
    {
      for (int i=0; i<_ncomp; i++)
        b[i] += __U[i].m() * N.m();
    }
};

/**********************************************************/
/**********************************************************/

MeshInterpolator::MeshInterpolator() : _MA(NULL), _DI(NULL)
{
}

/**********************************************************/

MeshInterpolator::~MeshInterpolator()
{
  if (_DI)
  {
    delete _DI;
    _DI = NULL;
  }
  if (_MA)
  {
    delete _MA;
    _MA = NULL;
  }
}

/**********************************************************/

void MeshInterpolator::CheckCell(int oldNumber, int newNumber)
{
  assert(_Old->sleep(oldNumber) && _New->sleep(newNumber));

  int oc0 = _Old->child(oldNumber,0), nc0 = _New->child(newNumber,0);
  if (_Old->sleep(oc0) && _New->sleep(nc0))
  {
    for (int i=0; i<_Old->nchilds(oldNumber); i++)
    {
      int oi = _Old->child(oldNumber,i);
      int ni = _New->child(newNumber,i);
      _NewCellNumber[oi] = ni;
      CheckCell(oi,ni);
    }
  }
  else
  {
    for (int j=0; j<_Old->nchilds(oldNumber); j++)
    {
      int ocj = _Old->child(oldNumber,j);
      int ncj = _New->child(newNumber,j);
      _NewCellNumber[ocj] = ncj;
      for (int i=0; i<_Old->nodes_per_cell(ocj); i++)
      {
        int oi = _Old->vertex_of_cell(ocj,i);
        _NewNodeNumber[oi] = _New->vertex_of_cell(ncj,i);
      }
    }
    if (!_VecNew.empty() && _Old->sleep(oc0) && !_New->sleep(nc0))
    {
      _ToBeRefNew.insert(newNumber);
    }
    else if (!_VecOld.empty() && !_Old->sleep(oc0) && _New->sleep(nc0))
    {
      _ToBeRef.insert(oldNumber);
    }
  }
}

/**********************************************************/

void MeshInterpolator::Coarsen(int newNumber)
{
  for (int i=0; i<_New->nchilds(newNumber); i++)
  {
    int ci = _New->child(newNumber,i);
    assert(_New->sleep(ci));
    if (_New->sleep(_New->child(ci,0)))
    {
      Coarsen(ci);
    }
  }
  int dim = _New->dimension();
  for (int c=0; c<_New->nchilds(newNumber); c++)
  {
    int cell = _New->child(newNumber,c);
    int npc = _New->nodes_per_cell(cell);
    for (int s=0; s<_VecInt.size(); s++)
    {
      int order = _VecInt[s].second;
      if (order==0)
      {
        _VecInt[s].first.zero_node(cell);
        double w = 1.;
        if (_average)
        {
          w /= static_cast<double>(npc);
        }
        for (int i=0; i<npc; i++)
        {
          int ci = _New->child(cell,i);
          _VecInt[s].first.add_node(cell,w,ci);
        }
      }
      else
      {
        if (order==1)
        {
          for (int i=0; i<npc; i++)
          {
            int ci = _New->child(cell,i);
            int k = _New->vertex_of_cell(cell,i);
            for (int j=0; j<npc; j++)
            {
              int l = _New->vertex_of_cell(ci,j);
              double w = _wq1(i,j);
              _VecInt[s].first.add_node(k,w,l);
            }
          }
        }
        else if (order==2)
        {
          int nind = (dim==2) ? 9 : 27;
          for (int i=0; i<nind; i++)
          {
            pair<int,int> indexI = _iq2[i];
            int ci = _New->child(newNumber,indexI.first);
            int k = _New->vertex_of_cell(ci,indexI.second);
            for (int j=0; j<nind; j++)
            {
              pair<int,int> indexJ = _iq2[j];
              int cj = _New->child(cell,indexJ.first);
              int l = _New->vertex_of_cell(cj,indexJ.second);
              double w = _wq2[c](i,j);
              _VecInt[s].first.add_node(k,w,l);
            }
          }
        }
        for (int i=0; i<npc; i++)
        {
          int ci = _New->child(cell,i);
          for (int j=0; j<npc; j++)
          {
            if (j!=i)
            {
              int l = _New->vertex_of_cell(ci,j);
              _VecInt[s].first.zero_node(l);
            }
          }
        }
      }
    }
  }
}

/**********************************************************/

void MeshInterpolator::Distribute(int oldNumber, int newNumber)
{
  assert(_Old->sleep(oldNumber) && _New->sleep(newNumber));

  int oc0 = _Old->child(oldNumber,0), nc0 = _New->child(newNumber,0);
  if (_Old->sleep(oc0) && _New->sleep(nc0))
  {
    for (int i=0; i<_Old->nchilds(oldNumber); i++)
    {
      int oi = _Old->child(oldNumber,i);
      int ni = _New->child(newNumber,i);
      Distribute(oi,ni);
    }
  }
  else if (!_Old->sleep(oc0) && _New->sleep(nc0))
  {
    Coarsen(newNumber);
  }
  else if (_Old->sleep(oc0) || _New->sleep(nc0))
  {
    cerr << "Das darf gar nicht passieren!!!" << endl;
    abort();
  }
}

/**********************************************************/

void MeshInterpolator::InitIndizes(int dim)
{
  int sizeq2 = static_cast<int>(pow(3.,static_cast<double>(dim)));
  _iq2.resize(sizeq2);
  if (dim==2)
  {
    _iq2[0] = make_pair(0,0);
    _iq2[1] = make_pair(0,1);
    _iq2[2] = make_pair(1,1);
    _iq2[3] = make_pair(0,3);
    _iq2[4] = make_pair(0,2);
    _iq2[5] = make_pair(1,2);
    _iq2[6] = make_pair(3,3);
    _iq2[7] = make_pair(2,3);
    _iq2[8] = make_pair(2,2);
  }
  else if (dim==3)
  {
    _iq2[ 0] = make_pair(0,0);
    _iq2[ 1] = make_pair(0,1);
    _iq2[ 2] = make_pair(1,1);
    _iq2[ 3] = make_pair(0,3);
    _iq2[ 4] = make_pair(0,2);
    _iq2[ 5] = make_pair(1,2);
    _iq2[ 6] = make_pair(3,3);
    _iq2[ 7] = make_pair(2,3);
    _iq2[ 8] = make_pair(2,2);
    _iq2[ 9] = make_pair(0,4);
    _iq2[10] = make_pair(0,5);
    _iq2[11] = make_pair(1,5);
    _iq2[12] = make_pair(0,7);
    _iq2[13] = make_pair(0,6);
    _iq2[14] = make_pair(1,6);
    _iq2[15] = make_pair(3,7);
    _iq2[16] = make_pair(2,7);
    _iq2[17] = make_pair(2,6);
    _iq2[18] = make_pair(4,4);
    _iq2[19] = make_pair(4,5);
    _iq2[20] = make_pair(5,5);
    _iq2[21] = make_pair(4,7);
    _iq2[22] = make_pair(4,6);
    _iq2[23] = make_pair(5,6);
    _iq2[24] = make_pair(7,7);
    _iq2[25] = make_pair(6,7);
    _iq2[26] = make_pair(6,6);
  }
}

/**********************************************************/

void MeshInterpolator::InitInterpolationWeights(int dim)
{
  // 1D-Gewichte fuer die automatische Berechnung der Q2-Gewichte
  nmatrix<double> w1d(3,5,0.);
  w1d(0,0) = 1.    ; w1d(0,1) = 0.375; w1d(0,3) = -0.125;
  w1d(1,1) = 0.75  ; w1d(1,2) = 1.   ; w1d(1,3) = 0.75  ;
  w1d(2,1) = -0.125; w1d(2,3) = 0.375; w1d(2,4) = 1.    ;
  
  int sizeq1 = static_cast<int>(pow(2.,static_cast<double>(dim)));
  int sizeq2 = static_cast<int>(pow(3.,static_cast<double>(dim)));
  _wq1.resize(sizeq1,sizeq1);
  _wq1.zero();

  _wq2.resize(sizeq1);
  for (int i=0; i<sizeq1; i++)
  {
    _wq2[i].resize(sizeq2,sizeq2);
    _wq2[i].zero();
  }

  if (dim==2)
  {
    // Q1-Gewichte
    _wq1(0,1) = 0.5  ; _wq1(0,2) = 0.25 ; _wq1(0,3) = 0.5  ;
    _wq1(1,0) = 0.5  ; _wq1(1,2) = 0.5  ; _wq1(1,3) = 0.25 ;
    _wq1(2,0) = 0.25 ; _wq1(2,1) = 0.5  ; _wq1(2,3) = 0.5  ;
    _wq1(3,0) = 0.5  ; _wq1(3,1) = 0.25 ; _wq1(3,2) = 0.5  ;

    // Q2-Gewichte
    for (int l=0; l<sizeq2; l++)
    {
      int lx = l%3;
      int ly = l/3;
      for (int i=0; i<3; i++)
      {
        for (int j=0; j<3; j++)
        {
          int m = i*3+j;
          _wq2[0](l,m) = w1d(ly,i)   * w1d(lx,j);
          _wq2[1](l,m) = w1d(ly,i)   * w1d(lx,j+2);
          _wq2[2](l,m) = w1d(ly,i+2) * w1d(lx,j+2);
          _wq2[3](l,m) = w1d(ly,i+2) * w1d(lx,j);
          for (int c=0; c<4; c++)
          {
            if (_wq2[c](l,m)==1.)
            {
              _wq2[c](l,m)=0.;
            }
          }
        }
      }
    }
  }
  else if (dim==3)
  {
    // Q1-Gewichte
    _wq1(0,1) = 0.5   ; _wq1(0,2) = 0.25  ; _wq1(0,3) = 0.5   ; _wq1(0,4) = 0.5   ; _wq1(0,5) = 0.25  ; _wq1(0,6) = 0.125 ; _wq1(0,7) = 0.25  ;
    _wq1(1,0) = 0.5   ; _wq1(1,2) = 0.5   ; _wq1(1,3) = 0.25  ; _wq1(1,4) = 0.25  ; _wq1(1,5) = 0.5   ; _wq1(1,6) = 0.25  ; _wq1(1,7) = 0.125 ;
    _wq1(2,0) = 0.25  ; _wq1(2,1) = 0.5   ; _wq1(2,3) = 0.5   ; _wq1(2,4) = 0.125 ; _wq1(2,5) = 0.25  ; _wq1(2,6) = 0.5   ; _wq1(2,7) = 0.25  ;
    _wq1(3,0) = 0.5   ; _wq1(3,1) = 0.25  ; _wq1(3,2) = 0.5   ; _wq1(3,4) = 0.25  ; _wq1(3,5) = 0.125 ; _wq1(3,6) = 0.25  ; _wq1(3,7) = 0.5   ;
    _wq1(4,0) = 0.5   ; _wq1(4,1) = 0.25  ; _wq1(4,2) = 0.125 ; _wq1(4,3) = 0.25  ; _wq1(4,5) = 0.5   ; _wq1(4,6) = 0.25  ; _wq1(4,7) = 0.5   ;
    _wq1(5,0) = 0.25  ; _wq1(5,1) = 0.5   ; _wq1(5,2) = 0.25  ; _wq1(5,3) = 0.125 ; _wq1(5,4) = 0.5   ; _wq1(5,6) = 0.5   ; _wq1(5,7) = 0.25  ;
    _wq1(6,0) = 0.125 ; _wq1(6,1) = 0.25  ; _wq1(6,2) = 0.5   ; _wq1(6,3) = 0.25  ; _wq1(6,4) = 0.25  ; _wq1(6,5) = 0.5   ; _wq1(6,7) = 0.5   ;
    _wq1(7,0) = 0.25  ; _wq1(7,1) = 0.125 ; _wq1(7,2) = 0.25  ; _wq1(7,3) = 0.5   ; _wq1(7,4) = 0.5   ; _wq1(7,5) = 0.25  ; _wq1(7,6) = 0.5   ;

    // Q2-Gewichte
    for (int l=0; l<sizeq2; l++)
    {
      int lx = l%3;
      int ly = (l%9)/3;
      int lz = l/9;
      for (int i=0; i<3; i++)
      {
        for (int j=0; j<3; j++)
        {
          for (int k=0; k<3; k++)
          {
            int m = i*9+j*3+k;
            _wq2[0](l,m) = w1d(lz,i)   * w1d(ly,j)   * w1d(lx,k);
            _wq2[1](l,m) = w1d(lz,i)   * w1d(ly,j)   * w1d(lx,k+2);
            _wq2[2](l,m) = w1d(lz,i)   * w1d(ly,j+2) * w1d(lx,k+2);
            _wq2[3](l,m) = w1d(lz,i)   * w1d(ly,j+2) * w1d(lx,k);
            _wq2[4](l,m) = w1d(lz,i+2) * w1d(ly,j)   * w1d(lx,k);
            _wq2[5](l,m) = w1d(lz,i+2) * w1d(ly,j)   * w1d(lx,k+2);
            _wq2[6](l,m) = w1d(lz,i+2) * w1d(ly,j+2) * w1d(lx,k+2);
            _wq2[7](l,m) = w1d(lz,i+2) * w1d(ly,j+2) * w1d(lx,k);
            for (int c=0; c<8; c++)
            {
              if (_wq2[c](l,m)==1.)
              {
                _wq2[c](l,m) = 0.;
              }
            }
          }
        }
      }
    }
  }
}

/**********************************************************/

void MeshInterpolator::RefineAndInterpolate(HierarchicalMesh* Mesh,
    vector<pair<GlobalVector,int> >& u, const IntSet& refine,
    vector<vector<bool> >& done)
{
  IntVector coarse(0), childs(0);
  int oldcells = Mesh->ncells();
  for (IntSet::const_iterator p=refine.begin(); p!=refine.end(); p++)
  {
    for (int j=0; j<Mesh->nchilds(*p); j++)
    {
      childs.push_back(Mesh->child(*p,j));
    }
  }
  Mesh->refine(childs,coarse);

  int nn = Mesh->nnodes(), nc = Mesh->ncells();
  for (int s=0; s<u.size(); s++)
  {
    if(u[s].second>0)
    {
      done[s].resize(nn,false);
      u[s].first.resize(nn,0.);
    }
    else
    {
      done[s].resize(nc,false);
      u[s].first.resize(nc,0.);
    }
  }
  IntSet fathers;
  for (int cell=oldcells; cell<Mesh->ncells(); cell++)
  {
    fathers.insert(Mesh->Vater(cell));
  }
  IntVector refined;
  for (IntSet::const_iterator p = fathers.begin(); p!=fathers.end(); p++)
  {
    refined.push_back(*p);
  }
  // nur ein Sicherheits-Check
  sort(childs.begin(),childs.end());
  sort(refined.begin(),refined.end());
  assert(childs==refined);

  int dim = Mesh->dimension();
  
  for (IntSet::const_iterator p=refine.begin(); p!=refine.end(); p++)
  {
    for (int s=0; s<u.size(); s++)
    {
      int order = u[s].second;
      for (int c=0; c<Mesh->nchilds(*p); c++)
      {
        int cell = Mesh->child(*p,c);
        int npc = Mesh->nodes_per_cell(cell);
        if (order==0)
        {
          for (int i=0; i<npc; i++)
          {
            int ci = Mesh->child(cell,i);
            u[s].first.add_node(ci,1.,cell);
            done[s][ci] = true;
          }
        }
        else
        {
          if (order==1)
          {
            for (int i=0; i<npc; i++)
            {
              int ci = Mesh->child(cell,i);
              int l = Mesh->vertex_of_cell(ci,i);
              for (int j=0; j<npc; j++)
              {
                int k = Mesh->vertex_of_cell(ci,j);
                if (!done[s][k])
                {
                  double w = _wq1(i,j);
                  u[s].first.add_node(k,w,l);
                }
              }
            }
          }
          else if (order==2)
          {
            int nind = (dim==2) ? 9 : 27;
            for (int i=0; i<nind; i++)
            {
              pair<int,int> indexI = _iq2[i];
              int ci = Mesh->child(*p,indexI.first);
              int l = Mesh->vertex_of_cell(ci,indexI.second);
              for (int j=0; j<nind; j++)
              {
                pair<int,int> indexJ = _iq2[j];
                int cj = Mesh->child(cell,indexJ.first);
                int k = Mesh->vertex_of_cell(cj,indexJ.second);
                if (!done[s][k])
                {
                  double w = _wq2[c](i,j);
                  u[s].first.add_node(k,w,l);
                }
              }
            }
          }
          for (int i=0; i<npc; i++)
          {
            int ci = Mesh->child(cell,i);
            for (int j=0; j<npc; j++)
            {
              if (j!=i)
              {
                int k = Mesh->vertex_of_cell(ci,j);
                done[s][k] = true;
              }
            }
          }
        }
      }
    }
  }
}

/**********************************************************/

void MeshInterpolator::AddVectorIntermediate(const GlobalVector& u, int order)
{
  _VecInt.push_back(make_pair(u,order));
}

/**********************************************************/

void MeshInterpolator::AddVectorOld(const GlobalVector& u, int order)
{
  _VecOld.push_back(make_pair(u,order));
  GetOriginalDiscretization()->HNAverage(_VecOld.back().first);
}

/**********************************************************/

void MeshInterpolator::AddVectorNew(const GlobalVector& u, int order)
{
  _VecNew.push_back(make_pair(u,order));
  GetDiscretization()->HNAverage(_VecNew.back().first);
}

/**********************************************************/

void MeshInterpolator::AddCellVectorOld(const GlobalVector& u)
{
  GlobalVector cu(u.ncomp(),_Old->ncells());
  const IntVector& celll2g = GetOriginalMeshAgent()->Celll2g();
  for(int i=0; i<u.n(); i++)
  {
    for(int c=0; c<u.ncomp(); c++)
    {
      cu(celll2g[i],c) = u(i,c);
    }
  }
  _VecOld.push_back(make_pair(cu,0));
}

/**********************************************************/

void MeshInterpolator::AddCellVectorNew(const GlobalVector& u)
{
  GlobalVector cu(u.ncomp(),_New->ncells());
  const IntVector& celll2g = GetMeshAgent()->Celll2g();
  for(int i=0; i<u.n(); i++)
  {
    for(int c=0; c<u.ncomp(); c++)
    {
      cu(celll2g[i],c) = u(i,c);
    }
  }
  _VecNew.push_back(make_pair(cu,0));
}

/**********************************************************/

void MeshInterpolator::BasicInit(DiscretizationInterface* DI, MeshAgentInterface* MA, const string& name)
{
  _name = name;

  // Original-Solver und -MeshAgent speichern
  MeshAgent *OMA = dynamic_cast<MeshAgent *>(MA);
  assert(OMA);
  _OMA = OMA;

  int dim = GetOriginalMeshAgent()->GetMesh(0)->dimension();
  _ODI = DI;

  // neuen MeshAgent anlegen
  _MA = new MeshAgent;

  GetMeshAgent()->GetShapes2d() = GetOriginalMeshAgent()->GetShapes2d();
  GetMeshAgent()->GetShapes3d() = GetOriginalMeshAgent()->GetShapes3d();
  GetMeshAgent()->BasicInit(_name+".gup",dim,0,0);
    
  // neue Discretization anlegen
  const Q2* Q2DP = dynamic_cast<const Q2*>(GetOriginalDiscretization());
  if (Q2DP)
  {
    if (dim==2)
    {
      _DI = new Q22d;
    }
    else
    {
      _DI = new Q23d;
    }
  }
  else
  {
    if (dim==2)
    {
      _DI = new Q12d;
    }
    else
    {
      _DI = new Q13d;
    }
  }
  GetDiscretization()->BasicInit(NULL);
  InitIndizes(dim);
  InitInterpolationWeights(dim);
  ReInit();
}

/**********************************************************/

void Gascoigne::MeshInterpolator::ReInit()
{
  // Klassenvariablen initialisieren
  _BaseCells.clear();
  _ToBeRef.clear();
  _ToBeRefNew.clear();
  _NewNodeNumber.clear();
  _NewCellNumber.clear();
  _VecInt.clear();
  _VecOld.clear();
  _VecNew.clear();
  _average = false;

  GetMeshAgent()->read_gup(_name);

  GetDiscretization()->ReInit(GetMeshAgent()->GetMesh(0));

  // Pointer auf die HierarchicalMeshs holen
  _Old = GetOriginalMeshAgent()->GetHierarchicalMesh();
  _New = GetMeshAgent()->GetHierarchicalMesh();
  assert(_Old);
  assert(_New);
}

/**********************************************************/

void Gascoigne::MeshInterpolator::RefineNodeVector(GlobalVector& uNew, const GlobalVector& uOld)
{
  const Q2* DI = dynamic_cast<const Q2*>(GetOriginalDiscretization());
  if (DI)
  {
    AddVectorNew(uOld,2);
  }
  else
  {
    AddVectorNew(uOld,1);
  }

  vector<vector<bool> > doneOld(_VecOld.size(),vector<bool>(_Old->nnodes(),true)),doneNew(_VecNew.size(),vector<bool>(_New->nnodes(),true));
  HierarchicalMesh* Mesh;
  if (_Old->ncells()<_New->ncells())
  {
    cerr << "Only possible if new mesh is finer than old mesh" << endl;
  }
  Mesh = _Old;

  for (int c=0; c<Mesh->ncells(); c++)
  {
    if (Mesh->level(c)==0)
    {
      _BaseCells.insert(c);
    }
  }

  do
  {
    _NewNodeNumber.resize(_Old->nnodes(),-1);
    _NewCellNumber.resize(_Old->ncells(),-1);
    _ToBeRef.clear();
    _ToBeRefNew.clear();
    for (IntSet::const_iterator pbc = _BaseCells.begin(); pbc!=_BaseCells.end(); pbc++)
    {
      CheckCell(*pbc,*pbc);
    }
    if (!_ToBeRef.empty())
    {
      RefineAndInterpolate(_Old,_VecOld,_ToBeRef,doneOld);
    }
    if (!_ToBeRefNew.empty())
    {
      RefineAndInterpolate(_New,_VecNew,_ToBeRefNew,doneNew);
    }
  }
  while (!_ToBeRef.empty() || !_ToBeRefNew.empty());

  GetMeshAgent()->ClearCl2g();
  GetMeshAgent()->global_refine(0);
  GetDiscretization()->ReInit(GetMeshAgent()->GetMesh(0));

  uNew.zero();
  assert(_VecNew.size()==1);
  for (int i=0; i<uNew.n(); i++)
  {
    uNew.equ_node(i,_NewNodeNumber[i],_VecNew[0].first);
  }

  _Old = NULL;
  _New = NULL;
}

/**********************************************************/

void MeshInterpolator::InterpolateCellVector(GlobalVector& out, const GlobalVector& in)
{
  AddCellVectorNew(in);

  vector<vector<bool> > doneOld(_VecOld.size(),vector<bool>(_Old->ncells(),true)),doneNew(_VecNew.size(),vector<bool>(_New->ncells(),true));
  HierarchicalMesh* Mesh;
  if (_Old->ncells()<_New->ncells())
  {
    Mesh = _Old;
  }
  else
  {
    Mesh = _New;
  }
  _NewCellNumber.resize(_Old->ncells(),-1);
  for (int c=0; c<Mesh->ncells(); c++)
  {
    if (Mesh->level(c)==0)
    {
      _BaseCells.insert(c);
      _NewCellNumber[c] = c;
    }
  }

  do
  {
    _NewNodeNumber.resize(_Old->nnodes(),-1);
    _NewCellNumber.resize(_Old->ncells(),-1);
    _ToBeRef.clear();
    _ToBeRefNew.clear();
    for (IntSet::const_iterator pbc = _BaseCells.begin(); pbc!=_BaseCells.end(); pbc++)
    {
      CheckCell(*pbc,*pbc);
    }
    if (!_ToBeRef.empty())
    {
      RefineAndInterpolate(_Old,_VecOld,_ToBeRef,doneOld);
    }
    if (!_ToBeRefNew.empty())
    {
      RefineAndInterpolate(_New,_VecNew,_ToBeRefNew,doneNew);
    }
  }
  while (!_ToBeRef.empty() || !_ToBeRefNew.empty());

  GetMeshAgent()->ClearCl2g();
  GetMeshAgent()->global_refine(0);

  assert(_VecNew.size()==1);
  _VecInt.push_back(_VecNew[0]);

  _average = true;
  for (IntSet::const_iterator pbc = _BaseCells.begin(); pbc!=_BaseCells.end(); pbc++)
  {
    Distribute(*pbc,*pbc);
  }
  _average = false;

  assert(_VecInt.size()==1);
  out.ReInit(in.ncomp(),GetOriginalMeshAgent()->GetMesh(0)->ncells());
  out.zero();

  const IntVector& cl2g = GetOriginalMeshAgent()->Celll2g();
  for (int i=0; i<out.n(); i++)
  {
    int n = cl2g[i];
    out.equ_node(i,_NewCellNumber[n],_VecInt[0].first);
  }

  _Old = NULL;
  _New = NULL;
}

/**********************************************************/

void MeshInterpolator::RhsForProjection(GlobalVector& f, const GlobalVector& u)
{
  const Q2* DI = dynamic_cast<const Q2*>(GetOriginalDiscretization());
  if (DI)
  {
    AddVectorNew(u,2);
  }
  else
  {
    AddVectorNew(u,1);
  }

  vector<vector<bool> > doneOld(_VecOld.size(),vector<bool>(_Old->nnodes(),true)),doneNew(_VecNew.size(),vector<bool>(_New->nnodes(),true));
  HierarchicalMesh* Mesh;
  if (_Old->ncells()<_New->ncells())
  {
    Mesh = _Old;
  }
  else
  {
    Mesh = _New;
  }
  for (int c=0; c<Mesh->ncells(); c++)
  {
    if (Mesh->level(c)==0)
    {
      _BaseCells.insert(c);
    }
  }

  do
  {
    _NewNodeNumber.resize(_Old->nnodes(),-1);
    _NewCellNumber.resize(_Old->ncells(),-1);
    _ToBeRef.clear();
    _ToBeRefNew.clear();
    for (IntSet::const_iterator pbc = _BaseCells.begin(); pbc!=_BaseCells.end(); pbc++)
    {
      CheckCell(*pbc,*pbc);
    }
    if (!_ToBeRef.empty())
    {
      RefineAndInterpolate(_Old,_VecOld,_ToBeRef,doneOld);
    }
    if (!_ToBeRefNew.empty())
    {
      RefineAndInterpolate(_New,_VecNew,_ToBeRefNew,doneNew);
    }
  }
  while (!_ToBeRef.empty() || !_ToBeRefNew.empty());

  GetMeshAgent()->ClearCl2g();
  GetMeshAgent()->global_refine(0);
  GetDiscretization()->ReInit(GetMeshAgent()->GetMesh(0));

  ProjectionRightHandSide PRHS(u.ncomp());
  
  GlobalVector _help(u.ncomp(),_New->nnodes(),0.);
  assert(_VecNew.size()==1);
  GetDiscretization()->AddNodeVector("U",&_VecNew[0].first);
  GetDiscretization()->Rhs(_help,PRHS,1.);
  GetDiscretization()->DeleteNodeVector("U");

  if (DI)
  {
    AddVectorIntermediate(_help,2);
  }
  else
  {
    AddVectorIntermediate(_help,1);
  }
  for (IntSet::const_iterator pbc = _BaseCells.begin(); pbc!=_BaseCells.end(); pbc++)
  {
    Distribute(*pbc,*pbc);
  }

  f.zero();
  assert(_VecInt.size()==1);
  for (int i=0; i<f.n(); i++)
  {
    f.equ_node(i,_NewNodeNumber[i],_VecInt[0].first);
  }
  GetOriginalDiscretization()->HNDistribute(f);

  _Old = NULL;
  _New = NULL;
}

/**********************************************************/
}

