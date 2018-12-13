/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the Gascoigne 3D authors
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


#include  <list>
#include  <iomanip> 

#include  "giota.h"
#include  "stdsolver.h"

#include  "pointmatrix.h"
#include  "pointilu.h"

#include  "dynamicblockmatrix.h"
#include  "dynamicblockilu.h"

/*--------------------------------*/
#ifdef __WITH_THREADS__
#include <omp.h>
#include  "threadilu.h"
#endif
/*--------------------------------*/

#include  "sparseblockilu.h"
#include  "fmatrixblock.h"
#include  "cfdblock3d.h"

/*--------------------------------*/
#ifdef __WITH_UMFPACK__
#include  "umfilu.h"
#endif
/*--------------------------------*/

#include  "ilupermutate.h"
#include  "cuthillmckee.h"
#include  "stopwatch.h"
#include  "gascoignevisualization.h"
#include  "backup.h"
#include  "visu_eps.h"
#include  "pi.h"

#include  "diracrighthandside.h"

#include  "q12d.h"
#include  "q22d.h"
#include  "q1gls2d.h"
#include  "q2gls2d.h"
#include  "q1lps2d.h"
#include  "q2lps2d.h"
#include  "q13d.h"
#include  "q23d.h"
#include  "q1gls3d.h"
#include  "q1lps3d.h"
#include  "q2lps3d.h"

#include  "faceq1.h"
#include  "faceq2.h"

#include "q2lps2dwithsecond.h"
#include "q22dwithsecond.h"
#include "q2lps3dwithsecond.h"
#include "q23dwithsecond.h"

#include  "glsequation.h"
#include  "lpsequation.h"

using namespace std;

/*-----------------------------------------*/
#ifdef __WITH_THREADS__
extern "C" void METIS_PartGraphRecursive(int *,int *,int *,int *,int *,int *,
					 int *,int *,int *,int *,int *);
extern "C" void METIS_PartGraphKway(int *,int *,int *,int *,int *,int *,
				    int *,int *,int *,int *,int *);
#endif
/*-----------------------------------------*/


namespace Gascoigne
{

StdSolver::StdSolver() : 
  _MP(NULL), _HM(NULL), _MAP(NULL), _MIP(NULL), _ZP(NULL), _FZP(NULL), _PDX(NULL), _NI(NULL),
  _distribute(true), _ndirect(1000), _directsolver(0), _discname("Q1"), _facediscname("none"),
  _matrixtype("point_node"), _PrimalSolve(1), _paramfile(NULL), _useUMFPACK(true)
// , omega_domain(0.) 
{
}

/*-----------------------------------------*/

StdSolver::~StdSolver()
{
  if(_MAP) delete _MAP; _MAP=NULL;
  if(_MIP) delete _MIP; _MIP=NULL;
  if(_ZP)  delete _ZP;  _ZP=NULL;
  if(_FZP) delete _FZP; _FZP=NULL;
}

/*-------------------------------------------------------*/

void StdSolver::_check_consistency(const Equation* EQ,const DiscretizationInterface* DI) const
{
  string eq = DI->GetName();

  bool glseq = false, glsdi = false;

  if (dynamic_cast<const GlsEquation*>(EQ))
  {
    glseq = true;
  }
  if (eq=="Q1Gls2d" || eq=="Q2Gls2d" || eq=="Q1Gls3d" || eq=="Q2Gls3d")
  {
    glsdi = true;
  }
  if(glseq && !glsdi)
  {
    cerr << "Warning: Discretization \"" << eq << "\" doesn't go with type of given Equation!" << endl;
  }
  else if(!glseq && glsdi)
  {
    cerr << "Error: Discretization \"" << eq << "\" doesn't go with type of given Equation!" << endl;
    abort();
  }


  bool lpseq = false, lpsdi = false;

  if(dynamic_cast<const LpsEquation*>(EQ))
  {
    lpseq = true;
  }
  if(eq=="Q1Lps2d" || eq=="Q2Lps2d" || eq=="Q1Lps3d" || eq=="Q2Lps3d")
  {
    lpsdi = true;
  }

  if(lpseq && !lpsdi)
  {
    cerr << "Warning: Discretization \"" << eq << "\" doesn't go with type of given Equation!" << endl;
  }
  else if(!lpseq && lpsdi)
  {
    cerr << "Error: Discretization \"" << eq << "\" doesn't go with type of given Equation!" << endl;
    abort();
  }
}

/*-------------------------------------------------------*/

void StdSolver::MatrixZero() const
{
  GetMatrix()->zero();
}

/*-------------------------------------------------------*/

void StdSolver::OutputSettings() const
{
  cout << "==================================================" << endl;
  cout << "Solver:                   " << GetName() << endl;
  cout << "Discretization:           " << GetDiscretization()->GetName()  << endl;
  if (GetFaceDiscretization())  
    cout << "FaceDiscretization:       " << GetFaceDiscretization()->GetName()  << endl;
  GetProblemDescriptor()->OutputSettings(cout);
  cout << "==================================================" << endl;
}

/*-------------------------------------------------------*/

void StdSolver::RegisterMatrix()
{
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);

  RegisterMatrix(EQ->GetNcomp());
}

/*-------------------------------------------------------*/

void StdSolver::RegisterMatrix(int ncomp)
{

#ifdef __WITH_UMFPACK__
  if (_useUMFPACK && _MAP!=NULL)
  {
    SimpleMatrix* SM = dynamic_cast<SimpleMatrix*>(GetMatrix());
    if ((SM && !_directsolver && _matrixtype!="point_node") || (!SM && _directsolver))
    {
      delete _MAP;
      _MAP = NULL;
    }
  }

  if (_useUMFPACK && _MIP!=NULL)
  {
    UmfIlu* UM = dynamic_cast<UmfIlu*>(GetIlu());
    if ((UM && !_directsolver) || (!UM && _directsolver))
    {
      delete _MIP;
      _MIP = NULL;
    }
  }
#endif

#ifdef __WITH_THREADS__
  if (__with_thread_ilu && _MIP!=NULL)
  {
    ThreadIlu* TIlu = dynamic_cast<ThreadIlu*>(GetIlu());
    if (!TIlu && __n_threads > 1)
    {
      delete _MIP; _MIP = NULL;
    }
    if (TIlu && __n_threads == 1)
    {
      delete _MIP; _MIP = NULL;
    }
  }
#endif

  if (_MAP==NULL)
    GetMatrixPointer() = NewMatrix(ncomp, _matrixtype);
  
  if (_MIP==NULL)
    GetIluPointer   () = NewIlu   (ncomp, _matrixtype);
}

/*-------------------------------------------------------*/

void StdSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  _PDX = &PDX;
  assert(_PDX);
  
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();

  if (EQ) _check_consistency(EQ,GetDiscretization());
}

/*-------------------------------------------------------*/

void StdSolver::SetDiscretization(DiscretizationInterface& DI, bool init)
{
  if(init)
  {
    DI.ReInit(GetMesh());
    DI.SetDataContainer(GetDiscretization()->GetDataContainer());
  }
  
  GetDiscretizationPointer() = &DI;
}

/*-------------------------------------------------------*/

#ifdef __WITH_THREADS__
void StdSolver::ThreadPartitionMesh()
{
  if(!__with_thread_ilu)
  {
    __n_threads = 1;
    return;
  }
  assert(GetMesh());
  const GascoigneMesh* M = dynamic_cast<const GascoigneMesh*> (GetMesh());
  assert(M);
  assert(M->HasPatch());

  int n_max_threads = omp_get_max_threads();
  int n_per_thread = __min_patches_per_thread;
  int __n          = M->npatches();
  __n_threads = std::max(1,std::min(n_max_threads, M->npatches()/n_per_thread));
  
  //nur wenns sinn macht
  if (__n_threads>1)
  {
    vector<int> adj1;
    vector<int> adj;
    adj.push_back(0);
    
    vector<vector<int> > node2patch(M->nnodes());
    for (int p=0;p<M->npatches();++p)
    {
      const vector<int>& ciop = M->CoarseIndices(p);
      for (int i=0;i<ciop.size();++i)
	node2patch[ciop[i]].push_back(p);
    }
// adjazenzliste erstellen
    for (int p=0;p<M->npatches();++p)
    {
      set<int> neighbors;
      const vector<int>& ciop = M->CoarseIndices(p);
      for (int i=0;i<ciop.size();++i)
	for (int n=0;n<node2patch[ciop[i]].size();++n)
	{
	  int neighbor = node2patch[ciop[i]][n];
	  neighbors.insert(neighbor);
      }
      for (set<int>::const_iterator it = neighbors.begin();it!=neighbors.end();++it)
      if (*it!=p)
	adj1.push_back(*it);
      adj.push_back(adj1.size());
    }
    assert(adj.size()==M->npatches()+1);
    assert(adj1.size()==adj[adj.size()-1]);
    
    int wgtflag    = 0;
    int numflag    = 0;
    int options[5] = {0,0,0,0,0};
    int edgecut    = -1;
    
    vector<int> patch_partition(M->npatches());
    
    if (__n_threads<8)
      METIS_PartGraphRecursive(&__n,&adj[0],&adj1[0],NULL,NULL,
			       &wgtflag,&numflag,&__n_threads,
			       &options[0],&edgecut,&patch_partition[0]);
    else
      METIS_PartGraphKway(&__n,&adj[0],&adj1[0], NULL, NULL,
			  &wgtflag,&numflag,&__n_threads,
			  &options[0],&edgecut,&patch_partition[0]);
  
    __thread_domain2node.clear();
    __thread_domain2node.resize(__n_threads);
    
    vector<set<int>    > __thread_domain2node_set(__n_threads);
    for (int p=0;p<M->npatches();++p)
    {
      const vector<int>& iop = *(M->IndicesOfPatch(p));
      __thread_domain2node_set[patch_partition[p]].insert(iop.begin(),iop.end());
    }
    for (int p=0;p<__n_threads;++p)
      for (set<int>::const_iterator it = __thread_domain2node_set[p].begin();
	   it!=__thread_domain2node_set[p].end();++it)
	__thread_domain2node[p].push_back(*it);
    
//The invers  mapping
    __thread_node2domain.clear();
    __thread_node2domain.resize(M->nnodes());
    for(int d=0; d<__n_threads; d++)
    {
      for(int n=0; n <  __thread_domain2node[d].size(); n++)
      {
	__thread_node2domain[__thread_domain2node[d][n]].push_back(make_pair(d,n));
      }
    }
  }//end of the case __n_threads > 1 sonst tun wir nix 
}
//Endof threads
#endif

/*-------------------------------------------------------*/


void StdSolver::NewMesh(const MeshInterface* mp)
{
  _MP = mp;
  assert(_MP);

  if(_MP->nnodes()<_ndirect) 
    {
      _directsolver=1;
    }
  else
    {
      _directsolver=0;
    }
  GetDiscretization()->ReInit(_MP);
  if (GetFaceDiscretization())
    GetFaceDiscretization()->ReInit(_MP);

  // 
  #ifdef __WITH_THREADS__
  ThreadPartitionMesh();
  #endif
}

/*-----------------------------------------*/

void StdSolver::SetDefaultValues(string discname, string matrixtype, int ndirect)
{
  _discname   = discname;
  _matrixtype = matrixtype;
  _ndirect    = ndirect;
}

/*-------------------------------------------------------*/

void StdSolver::BasicInit(const ParamFile* paramfile, const int dimension,const NumericInterface* NI)
{
  _paramfile = paramfile;
  _NI = NI;

  string xxx;

  DataFormatHandler DFH;
  DFH.insert("matrixtype" , &_matrixtype);
  DFH.insert("ndirect"    , &_ndirect);
  DFH.insert("useUMFPACK", &_useUMFPACK);
  DFH.insert("discname",    &_discname);
  DFH.insert("facediscname",    &_facediscname,"none");
  DFH.insert("disc", &xxx, "void");

#ifdef __WITH_THREADS__ 
  int default_min_per_thread = 450;
  if(dimension==3)
    default_min_per_thread = 150;
  //Default is approximately 4000 nodes per Thread, the value is the 
  //minumum number of Patches for each thread, e.g in 2d 450*9 \approx 4000   
  DFH.insert("min_patches_per_thread",  &__min_patches_per_thread,default_min_per_thread); 
  DFH.insert("with_thread_ilu", &__with_thread_ilu, false);
#endif

  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Solver");
 
  if(xxx!="void")
    {
      cout << "Expression 'disc' in ParamFile not longer valid !" << endl;
      abort();
    }

#ifdef __WITH_THREADS__
  if(__with_thread_ilu && _matrixtype!="block")
  {
    cerr << "Thread Ilu (with_thread_ilu) not available for matrixtype " << _matrixtype << "!" << endl;
    cerr << "Select matrixtype block if you want to use it." << endl;
    abort();
  }
#endif

  if (GetDiscretizationPointer()==NULL)
    GetDiscretizationPointer()     = NewDiscretization    (dimension, _discname);
  GetFaceDiscretizationPointer() = NewFaceDiscretization(dimension, _facediscname);
  assert(_ZP);

  GetDiscretization()->BasicInit(_paramfile);
  if (GetFaceDiscretization())
    GetFaceDiscretization()->BasicInit(_paramfile);
  

  _Dat.BasicInit(_paramfile);
  _PF.SetComponents(_Dat.GetPfilter());
}

/*-------------------------------------------------------*/

DiscretizationInterface* StdSolver::NewDiscretization(int dimension, const string& discname)
{
  if(_NI)
  {
    return _NI->NewDiscretization();
  }
  else
  {
  }
  if (dimension==2)
  {
    if      (discname=="Q1")               return new Q12d;
    else if (discname=="Q2")               return new Q22d;
    else if (discname=="Q1Gls")            return new Q1Gls2d;
    else if (discname=="Q2Gls")            return new Q2Gls2d;
    else if (discname=="Q1Lps")            return new Q1Lps2d;
    else if (discname=="Q2Lps")            return new Q2Lps2d;
    else if (discname=="Q2WithSecond")     return new Q22dWithSecond;
    else if (discname=="Q2LpsWithSecond")  return new Q2Lps2dWithSecond;
    else 
    {         
      cerr << " Solver::NewDiscretization()\tunknown discname=" << discname << endl;
      abort();
    }
  }
  else if (dimension==3)
  {
    if      (discname=="Q1")               return new Q13d;
    else if (discname=="Q2")               return new Q23d;
    else if (discname=="Q1Gls")            return new Q1Gls3d;
    else if (discname=="Q1Lps")            return new Q1Lps3d;
    else if (discname=="Q2Lps")            return new Q2Lps3d;
    else if (discname=="Q2WithSecond")     return new Q23dWithSecond;
    else if (discname=="Q2LpsWithSecond")  return new Q2Lps3dWithSecond;
    else 
    {         
      cerr << " Solver::NewDiscretization()\tunknown discname=" << discname << endl;
      abort();
    }
  }
  else
  {
    cerr << " Solver::NewDiscretization()\tdimension must either 2 or 3" << endl;
    abort();

  }
}

/*-------------------------------------------------------------*/

FaceDiscretization* StdSolver::NewFaceDiscretization(int dimension, const string& facediscname)
{
  if (facediscname=="none") return NULL;

  if (dimension==2)
    {
      if      (facediscname=="FaceQ1") return new FaceQ1<2>;
      else if (facediscname=="FaceQ2") return new FaceQ2<2>;
      else 
	{
	  cerr << " Solver::NewDiscretization()\tunknown facediscname=" << facediscname << endl;
	      abort();
	}
    }
  else if (dimension==3)
    {
      if      (facediscname=="FaceQ1") return new FaceQ1<3>;
      else if (facediscname=="FaceQ2") return new FaceQ2<3>;
      else 
	{
	  cerr << "Solver::NewFaceDiscretization()\tunknown facediscname=" << facediscname << endl;
	  abort();
	}
    }
  else
    {
      cerr << "Solver::NewFaceDiscretization()\tdimension must either 2 or 3" << endl;
      abort();
    }
}

/*-------------------------------------------------------------*/

MatrixInterface* StdSolver::NewMatrix(int ncomp, const string& matrixtype) 
{
  if( _directsolver || matrixtype=="point_node")
  {
    return new PointMatrix(ncomp,"node");
  }
  else if (matrixtype=="block")
  {
    if      (ncomp==1)  return new SparseBlockMatrix<FMatrixBlock<1> >;
    else if (ncomp==2)  return new SparseBlockMatrix<FMatrixBlock<2> >;
    else if (ncomp==3)  return new SparseBlockMatrix<FMatrixBlock<3> >;
    else if (ncomp==4)  return new SparseBlockMatrix<FMatrixBlock<4> >;
    else
    {
      cerr << "No SparseBlockMatrix for " << ncomp << "components." << endl;
      abort();
    }
  }
  else if (matrixtype=="dynamic")
  {
    if      (ncomp==1)  return new DynamicBlockMatrix<FMatrixBlock<1> >;
    else if (ncomp==2)  return new DynamicBlockMatrix<FMatrixBlock<2> >;
    else if (ncomp==3)  return new DynamicBlockMatrix<FMatrixBlock<3> >;
    else if (ncomp==4)  return new DynamicBlockMatrix<FMatrixBlock<4> >;
    else
    {
      cerr << "No SparseBlockMatrix for " << ncomp << "components." << endl;
      abort();
    }
  }
  else if (matrixtype=="component")
  {
    return new PointMatrix(ncomp,"component");
  }
  else if (matrixtype=="cfd")
  {
    if (ncomp==4)
    {
      return new SparseBlockMatrix<CFDBlock3d>;
    }
    else
    {
      cerr << "No SparseBlockMatrix for " << ncomp << "components." << endl;
      abort();
    }
  }
  else
  {
    cerr << "No such matrix type \"" << matrixtype<< "\"." << endl;
    abort();
  }
}

/*-------------------------------------------------------------*/

IluInterface* StdSolver::NewIlu(int ncomp, const string& matrixtype) 
{ 

#ifdef __WITH_UMFPACK__
  if(_directsolver && _useUMFPACK)
  {
    return new UmfIlu(GetMatrix());
  }
#endif

  // analog zu NewMatrix muss hier auch _directsolver eingehen, 
  // sonst gibts aerger nachher beim 
  // GetIlu()->copy_entries(GetMatrix());
  if(_directsolver || matrixtype=="point_node")
  {
    return new PointIlu(ncomp,"node");
  }
  
  else if (matrixtype=="block")
  {

#ifdef __WITH_THREADS__
    if(__n_threads > 1)
    {
      return new ThreadIlu(ncomp);
    }
#endif 

    if      (ncomp==1)  return new SparseBlockIlu<FMatrixBlock<1> >;
    else if (ncomp==2)  return new SparseBlockIlu<FMatrixBlock<2> >;
    else if (ncomp==3)  return new SparseBlockIlu<FMatrixBlock<3> >;
    else if (ncomp==4)  return new SparseBlockIlu<FMatrixBlock<4> >;
    else
      {
        cerr << "No SparseBlockIlu for " << ncomp << "components." << endl;
        abort();
      }
  }
  else if (matrixtype=="dynamic")
   {
    if      (ncomp==1)  return new DynamicBlockIlu<FMatrixBlock<1> >;
    else if (ncomp==2)  return new DynamicBlockIlu<FMatrixBlock<2> >;
    else if (ncomp==3)  return new DynamicBlockIlu<FMatrixBlock<3> >;
    else if (ncomp==4)  return new DynamicBlockIlu<FMatrixBlock<4> >;
    else
      {
        cerr << "No DynamicBlockIlu for " << ncomp << "components." << endl;
        abort();
      }
  }
 else if (matrixtype=="component") 
  {
    return new PointIlu(ncomp,"component");
  }
  else if (matrixtype=="cfd")
  {
    if (ncomp==4)  
    {
      return new SparseBlockIlu<CFDBlock3d>;
    }
  }
  cerr << "No such matrix type \"" << matrixtype << "and ncomp \"." << ncomp << endl;
  abort();
}

/*-------------------------------------------------------*/

void StdSolver::ReInitMatrix() 
{
  GetDiscretization()->InitFilter(GetPfilter());
  SparseStructure SA;
  GetDiscretization()->Structure(&SA);
  
  if (GetFaceDiscretization())
    GetFaceDiscretization()->Structure(&SA);

  AddPeriodicNodes(&SA);

  GetMatrix()->ReInit(&SA);
  GetIlu()->ReInit(&SA);
}

/*-------------------------------------------------------*/

void StdSolver::AddPeriodicNodes(SparseStructure* SA)
{
  /*-------------------------------------------------------
  | Vereinigt die Nachbarschaften von zusammengehoerenden
  | Knoten auf den periodischen Raendern, so dass es
  | moeglich wird, die Kopplungen des einen Knoten zum
  | anderen zu schieben.
  -------------------------------------------------------*/
  
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntVector& iv_PeriodicColors = BM->GetPeriodicDataColors();

  const MeshInterface* p_mesh        = GetMesh();
  const GascoigneMesh* GMP           = dynamic_cast<const GascoigneMesh*>(p_mesh);
  assert(GMP);

  map<int,map<int,int> > mm_PeriodicPairs = GMP->GetBoundaryIndexHandler().GetPeriodicPairs();

  for(IntVector::const_iterator p_col   = iv_PeriodicColors.begin(); p_col!=iv_PeriodicColors.end(); p_col += 2)
  {
    int col = *p_col;

    IntSet is_neighbours1;
    IntSet is_neighbours2;

    for (map<int,int>::const_iterator p_pair = mm_PeriodicPairs[col].begin(); p_pair != mm_PeriodicPairs[col].end(); p_pair++)
    {
      // beide raender abgrasen und die kopplungen in columns1 und columns2 eintragen
      is_neighbours1 = SA->row(p_pair->first);
      is_neighbours2 = SA->row(p_pair->second);

      for (IntSet::const_iterator p_neighbour1 = is_neighbours1.begin(); p_neighbour1 != is_neighbours1.end(); p_neighbour1++)
      {
        SA->build_add(p_pair->second, *p_neighbour1  );
        SA->build_add(*p_neighbour1 , p_pair->second );
      }

      for (IntSet::const_iterator p_neighbour2 = is_neighbours2.begin(); p_neighbour2 != is_neighbours2.end(); p_neighbour2++)
      {
        SA->build_add(p_pair->first, *p_neighbour2 );
        SA->build_add(*p_neighbour2, p_pair->first );
      }
    }
  }
  SA->build_end();
}

/*-------------------------------------------------------*/

void StdSolver::RegisterVector(const VectorInterface& g) 
{
  _NGVA.Register(g);
}

/*-------------------------------------------------------*/

void StdSolver::ReInitVector(VectorInterface& dst)
{
  int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
  ReInitVector(dst,ncomp);
}

/*-------------------------------------------------------*/

void StdSolver::ReInitVector(VectorInterface& dst, int comp)
{
  int n = GetDiscretization()->n();
  int nc = GetDiscretization()->nc();

  // VectorInterface already registered ?
  //
  GhostVectorAgent::iterator p = _NGVA.find(dst);
  if (p==_NGVA.end()) 
    {
      _NGVA.Register(dst);
      p = _NGVA.find(dst);
    }
  assert(p!=_NGVA.end());
  
  // GlobalVector already registered ?
  //
  if (p->second==NULL) 
    {
      p->second = new GlobalVector;
    }

  // resize GlobalVector
  //
  p->second->ncomp()=comp;

  if(p->first.GetType()=="node")
    {
      p->second->reservesize(n);
    }
  else if(p->first.GetType()=="cell")
    {
      p->second->reservesize(nc);
    }
  else if(p->first.GetType()=="parameter")
    {
      p->second->reservesize(1);
    }
  else
  {
    cerr << "No such vector type: " << p->first.GetType() << endl;
    abort();
  }
}

/*-------------------------------------------------------*/

void StdSolver::Zero(VectorInterface& dst) const
{
  GetGV(dst).zero();
}

/*-----------------------------------------*/

double StdSolver::NewtonNorm(const VectorInterface& u) const
{
  return GetGV(u).norm_l8();
}

/*-----------------------------------------*/

void StdSolver::HNAverageData() const {
  GetDiscretization()->HNAverageData();
}
void StdSolver::HNZeroData() const {
  GetDiscretization()->HNZeroData();
}
void StdSolver::HNAverage(const VectorInterface& x) const {
  GetDiscretization()->HNAverage(const_cast<GlobalVector&>(GetGV(x)));
}
void StdSolver::HNZero(const VectorInterface& x) const {
  GetDiscretization()->HNZero(const_cast<GlobalVector&>(GetGV(x)));;
}
void StdSolver::HNDistribute(VectorInterface& x) const {
  if(GetDistribute())
  {
    GetDiscretization()->HNDistribute(GetGV(x));
  }
}

/*-------------------------------------------------------*/

void StdSolver::InterpolateSolution(VectorInterface& gu, const GlobalVector& uold) const
{
  GlobalVector& u = GetGV(gu);

  u.zero();
  GetDiscretization()->InterpolateSolution(u, uold);
  SubtractMean(gu);
}

/*-----------------------------------------*/

void StdSolver::residualgmres(VectorInterface& gy, const VectorInterface& gx, const VectorInterface& gb) const
{
  GlobalVector& y = GetGV(gy);
  const GlobalVector& b = GetGV(gb);

  vmulteq(gy,gx,1.);
  y.sadd(-1.,1.,b);
  SetBoundaryVectorZero(gy);
}

/*-----------------------------------------*/

void StdSolver::vmult(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  _vm.start();
  GetMatrix()->vmult(GetGV(gy),GetGV(gx),d);
  _vm.stop();
}

/*-----------------------------------------*/

void StdSolver::vmulteq(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  _vm.start();
  Zero(gy);
  _vm.stop();
  vmult(gy,gx,d);
}

/*-----------------------------------------*/

void StdSolver::MatrixResidual(VectorInterface& gy, const VectorInterface& gx, const VectorInterface& gb) const
{
  Equ(gy,1.,gb);
  vmult(gy,gx,-1.);
  SubtractMeanAlgebraic(gy);
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorZero(VectorInterface& gf) const
{
  GlobalVector& f = GetGV(gf);

  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntSet& Colors = BM->GetDirichletDataColors();
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM->GetDirichletDataComponents(col);
      GetDiscretization()->StrongDirichletVectorZero(f, col, comp);
    }
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVector(VectorInterface& gf) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const DirichletData*   DD = GetProblemDescriptor()->GetDirichletData();
  if(DD==NULL) 
  {
    if(BM->GetDirichletDataColors().size()!=0) 
    {
      cerr << "No DirichetData given but DirichetColors in ParamFile!" << endl;
      abort();
    }
    return;
  }
  SetBoundaryVectorStrong(gf,*BM,*DD);
}

/*-------------------------------------------------------*/

void StdSolver::SetPeriodicVector(VectorInterface& gf) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const PeriodicData*    PD = GetProblemDescriptor()->GetPeriodicData();

  if(PD==NULL)
  {
    if(BM->GetPeriodicDataColors().size()!=0)
    {
      cerr << "No PeriodicData given but PeriodicColors in ParamFile!" << endl;
      abort();
    }
    return;
  }

  SetPeriodicVectorStrong(gf,*BM,*PD);
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorStrong(VectorInterface& gf, const BoundaryManager& BM, const DirichletData& DD, double d) const
{
  GlobalVector& f = GetGV(gf);

  IntSet PrefCol = DD.preferred_colors();
  list<int> colors(BM.GetDirichletDataColors().begin(), 
		   BM.GetDirichletDataColors().end());
  
  for(IntSet::const_iterator p=PrefCol.begin();p!=PrefCol.end();p++)
    {
      int col = *p;
      colors.remove(col);
      colors.push_back(col);
    }
  for(list<int>::const_iterator p=colors.begin();p!=colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM.GetDirichletDataComponents(col);
      GetDiscretization()->StrongDirichletVector(f, DD, col, comp, d);
    }
}

/*-------------------------------------------------------*/

void StdSolver::SetPeriodicVectorStrong(VectorInterface& gf, const BoundaryManager& BM, const PeriodicData& PD, double d) const
{
  GlobalVector& f = GetGV(gf);

  list<int> periodic_cols(BM.GetPeriodicDataColors().begin(),
                            BM.GetPeriodicDataColors().end());
  for(list<int>::const_iterator p=periodic_cols.begin();p!=periodic_cols.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM.GetPeriodicDataComponents(col);
      GetDiscretization()->StrongPeriodicVector(f, PD, col, comp, d);
    }
}

/*-------------------------------------------------------*/

void StdSolver::SetPeriodicVectorZero(VectorInterface& gf) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntVector& iv_PeriodicColors = BM->GetPeriodicDataColors();

  GlobalVector& f = GetGV(gf);
  const MeshInterface* p_mesh        = GetMesh();
  const GascoigneMesh* GMP           = dynamic_cast<const GascoigneMesh*>(p_mesh);
  assert(GMP);

  map<int,map<int,int> > mm_PeriodicPairs = GMP->GetBoundaryIndexHandler().GetPeriodicPairs();

  for(IntVector::const_iterator p_col=iv_PeriodicColors.begin(); p_col!=iv_PeriodicColors.end(); )
  {
    int col  = *p_col++;
    *p_col++;

    const IntVector& iv_PeriodicComponents = BM->GetPeriodicDataComponents(col);

    for (map<int,int>::const_iterator p_pair = mm_PeriodicPairs[col].begin(); p_pair != mm_PeriodicPairs[col].end(); p_pair++)
    {
      for (IntVector::const_iterator p_comp=iv_PeriodicComponents.begin();p_comp!=iv_PeriodicComponents.end();p_comp++)
      {
        f(p_pair->second,*p_comp) = .5*(f(p_pair->second,*p_comp) + f(p_pair->first,*p_comp));
        f(p_pair->first ,*p_comp) = f(p_pair->second,*p_comp);
      }
    }
  }
}

/*-------------------------------------------------------*/

void StdSolver::smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const
{
  _il.start();
  double omega = GetSolverData().GetOmega();
  
  for(int iter=0; iter<niter; iter++)
    {
      if (GetSolverData().GetLinearSmooth()=="ilu")
	{
	  MatrixResidual(h,x,y);
	  GetIlu()->solve(GetGV(h));
	  Add(x,omega,h);
	}
      else if (GetSolverData().GetLinearSmooth()=="jacobi")
	{
	  MatrixResidual(h,x,y);
	  GetMatrix()->Jacobi(GetGV(h));
	  Add(x,omega,h);
	}
      else if (GetSolverData().GetLinearSmooth()=="richardson")
	{
	  MatrixResidual(h,x,y);
	  Add(x,omega,h);
	}
      else if (GetSolverData().GetLinearSmooth()=="none")
	{}
      else
	{
	  cerr << "Smoother: " << GetSolverData().GetLinearSmooth() << " not valid!\n";
	  abort();
	}
      SubtractMean(x);
    }
  _il.stop();
}

/*-------------------------------------------------------*/

void StdSolver::smooth_pre(VectorInterface& x, const VectorInterface& y, VectorInterface& help) const
{
  int niter = GetSolverData().GetIterPre();
  smooth(niter,x,y,help);
}

/*-------------------------------------------------------*/

void StdSolver::smooth_exact(VectorInterface& x, const VectorInterface& y, VectorInterface& help) const
{
#ifdef __WITH_UMFPACK__
  if(_directsolver&&_useUMFPACK)
    {
      _so.start();
      UmfIlu* UM = dynamic_cast<UmfIlu*>(GetIlu());
      assert(UM);
      UM->Solve(GetGV(x),GetGV(y));
      _so.stop();
    }
  else
#endif
    {
      int niter = GetSolverData().GetIterExact();
      smooth(niter,x,y,help);
    }
}

/*-------------------------------------------------------*/

void StdSolver::smooth_post(VectorInterface& x, const VectorInterface& y, VectorInterface& help) const
{
  int niter = GetSolverData().GetIterPost();
  smooth(niter,x,y,help);
}

/*-------------------------------------------------------*/

void StdSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  _re.start();

  HNAverage(gx);
  HNAverageData();

  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  GetDiscretization()->Form(GetGV(gy),GetGV(gx),*EQ,d);

  // Face
  if (GetFaceDiscretization())
    {
      const FaceEquation* FEQ = GetProblemDescriptor()->GetFaceEquation();
      assert(FEQ);
      GetFaceDiscretization()->FaceForm(GetGV(gy),GetGV(gx),*FEQ,d);
    }  
  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
  {
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetDiscretization()->BoundaryForm(GetGV(gy),GetGV(gx),BM->GetBoundaryEquationColors(),*BE,d);
  }
  HNZero(gx);
  HNZeroData();
  HNDistribute(gy);
  SubtractMeanAlgebraic(gy);

  _re.stop();
}

/*-------------------------------------------------------*/

void StdSolver::AdjointForm(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  GlobalVector&       y = GetGV(gy);
  const GlobalVector& x = GetGV(gx);
  
  HNAverage(gx);
  HNAverageData();

  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  GetDiscretization()->AdjointForm(y,x,*EQ,d);

  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
  {
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetDiscretization()->BoundaryForm(y,x,BM->GetBoundaryEquationColors(),*BE,d);
  }

  HNZero(gx);
  HNZeroData();
  HNDistribute(gy);
  SubtractMeanAlgebraic(gy);
}

/*-------------------------------------------------------*/

void StdSolver::BoundaryInit(VectorInterface& Gu)  const
{
  GlobalVector& u = GetGV(Gu);
  const DirichletData* DD = GetProblemDescriptor()->GetDirichletData();
  if(DD==NULL)
    {
      u.zero();
      cerr << "StdSolver::BoundaryInit():\t No DirichletData !!\n";
      abort();
    }
  
  int color = 0;
  int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
  DoubleVector y(ncomp);
  
  int dim = GetMesh()->dimension();
  for (int ind=0; ind<GetMesh()->nnodes(); ind++)
    {
      if (dim==2) (*DD)(y,GetMesh()->vertex2d(ind),color);
      else        (*DD)(y,GetMesh()->vertex3d(ind),color);
      
      for (int c=0; c<u.ncomp(); c++)
        {
          u(ind,c) = y[c];
        }
    }
}

/*-------------------------------------------------------*/

void StdSolver::SolutionInit(VectorInterface& Gu)  const
{
  GlobalVector& u = GetGV(Gu);
  const DomainInitialCondition* u0 = dynamic_cast<const DomainRightHandSide *>(GetProblemDescriptor()->GetInitialCondition());
  if(u0==NULL)
    {
      u.zero();
      return;
    }

  assert(u.ncomp()==GetProblemDescriptor()->GetEquation()->GetNcomp());
  
  for (int ind=0; ind<GetMesh()->nnodes(); ind++)
    {
      if (GetMesh()->dimension()==2)
        {
          for (int c=0; c<u.ncomp(); c++)
            {
              u(ind,c) = (*u0)(c,GetMesh()->vertex2d(ind));
            }
        }
            else
        {
          for (int c=0; c<u.ncomp(); c++)
            {
              u(ind,c) = (*u0)(c,GetMesh()->vertex3d(ind));
            }
        }
    }
}

/*-------------------------------------------------------*/

void StdSolver::ComputeError(const VectorInterface& u, GlobalVector& err) const
{
  if(GetProblemDescriptor()->GetExactSolution()==NULL) return;
  HNAverage(u);
  GetDiscretization()->ComputeError(GetGV(u),err,GetProblemDescriptor()->GetExactSolution());
  HNZero(u);
}

/*-------------------------------------------------------*/

void StdSolver::AssembleError(GlobalVector& eta, const VectorInterface& u, GlobalVector& err) const
{
  if(GetProblemDescriptor()->GetExactSolution()==NULL) return;
  HNAverage(u);
  GetDiscretization()->AssembleError(eta,GetGV(u),err,GetProblemDescriptor()->GetExactSolution());
  HNZero(u);
}

/*-------------------------------------------------------*/

double StdSolver::ComputeFunctional(VectorInterface& gf, const VectorInterface& gu, const Functional* FP)
{
  VectorInterface gh("XXX");
  ReInitVector(gh);
  double val = 0.;

  const DomainFunctional*   DFP = dynamic_cast<const DomainFunctional*>(FP);
  const BoundaryFunctional* BFP = dynamic_cast<const BoundaryFunctional*>(FP);
  const ResidualFunctional* RFP = dynamic_cast<const ResidualFunctional*>(FP);
  const PointFunctional*   NPFP = dynamic_cast<const PointFunctional*>(FP);

  if      (DFP) val = ComputeDomainFunctional(gf,gu,gh,DFP);
  else if (BFP) val = ComputeBoundaryFunctional(gf,gu,gh,BFP);
  else if(RFP)  val = ComputeResidualFunctional(gf,gu,gh,RFP);
  else if(NPFP) val = ComputePointFunctional(gf,gu,gh,NPFP);
  else
    {
      cerr << "Functional must be either of type DomainFunctional, BoundaryFunctional or PointFunctional!!!" << endl;
      abort();
    }
  DeleteVector(gh);
  return val;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeBoundaryFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, const BoundaryFunctional* FP) const
{
  HNAverage(gu);
  HNAverageData();
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  double J = GetDiscretization()->ComputeBoundaryFunctional(GetGV(gu),BM->GetBoundaryFunctionalColors(),*FP);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeDomainFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, const DomainFunctional* FP) const
{
  HNAverage(gu);
  HNAverageData();
  double J = GetDiscretization()->ComputeDomainFunctional(GetGV(gu),*FP);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputePointFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, const PointFunctional* FP) const
{
  HNAverage(gu);
  HNAverageData();
  double J = GetDiscretization()->ComputePointFunctional(GetGV(gu),*FP);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeResidualFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, const ResidualFunctional* FP) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();

  HNAverage(gu);
  HNAverageData();
  Zero(gf);
  Rhs(gf);


  Form(gf,gu,-1.);

  const DirichletData* ABD = FP->GetDirichletData();
  assert(ABD);

  Zero(gz);
  SetBoundaryVectorStrong(gz,*BM,*ABD);
  
  HNAverage(gz);

  double J = ScalarProduct(gz,gf);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

void StdSolver::EvaluateCellRightHandSide(VectorInterface& f, const DomainRightHandSide& CF, double d) const
{
  assert(f.GetType()=="cell");
  HNAverageData();
  
  GetDiscretization()->EvaluateCellRightHandSide(GetGV(f),CF,d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void StdSolver::EvaluateBoundaryCellRightHandSide(VectorInterface& f, const BoundaryRightHandSide& CF,
    const BoundaryManager& BM, double d) const
{
  assert(f.GetType()=="cell");
  HNAverageData();
  
  GetDiscretization()->EvaluateBoundaryCellRightHandSide(GetGV(f),BM.GetBoundaryRightHandSideColors(),CF,d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void StdSolver::EvaluateParameterRightHandSide(VectorInterface& f, const DomainRightHandSide& CF, double d) const
{
  assert(f.GetType()=="parameter");
  HNAverageData();
  
  GetDiscretization()->EvaluateParameterRightHandSide(GetGV(f),CF,d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void StdSolver::EvaluateBoundaryParameterRightHandSide(VectorInterface& f, const BoundaryRightHandSide& CF,
    const BoundaryManager& BM, double d) const
{
  assert(f.GetType()=="parameter");
  HNAverageData();
  
  GetDiscretization()->EvaluateBoundaryParameterRightHandSide(GetGV(f),BM.GetBoundaryRightHandSideColors(),CF,d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void StdSolver::InterpolateDomainFunction(VectorInterface&  f, const DomainFunction& DF) const
{
  HNAverageData();
  
  if(f.GetType()=="node")
  {
    GetDiscretization()->InterpolateDomainFunction(GetGV(f),DF);
  }
  else if(f.GetType()=="cell")
  {
    GetDiscretization()->InterpolateCellDomainFunction(GetGV(f),DF);
  }
  else
  {
    cerr << "No such vector type: " << f.GetType() << endl;
    abort();
  }

  HNZeroData();
}

/*-------------------------------------------------------*/

void StdSolver::Rhs(VectorInterface& gf, double d) const
{
  GlobalVector& f = GetGV(gf);
  HNAverageData();

  const Application* RHS  = GetProblemDescriptor()->GetRightHandSide();
  const BoundaryRightHandSide* NRHS = GetProblemDescriptor()->GetBoundaryRightHandSide();

  if(RHS)
    {
       bool done=false;
       const DomainRightHandSide *DRHS = dynamic_cast<const DomainRightHandSide *>(RHS);
       if(DRHS)
       {
         GetDiscretization()->Rhs(f,*DRHS,d);
         done = true;
       }
       const DiracRightHandSide *NDRHS = dynamic_cast<const DiracRightHandSide *>(RHS);
       if(NDRHS)
       {
         GetDiscretization()->DiracRhs(f,*NDRHS,d);
         done =true;
       }
       if(!done)
       {
         cerr << "RightHandSide should be either of type DomainRightHandSide or DiracRightHandSide!!!" << endl;
         abort();
       }
    }
  
  if(NRHS)
    {
      assert(NRHS->GetNcomp()==f.ncomp());
      const BoundaryManager*  BM   = GetProblemDescriptor()->GetBoundaryManager();
      GetDiscretization()->BoundaryRhs(f,BM->GetBoundaryRightHandSideColors(),*NRHS,d);	  
    }
    
  HNZeroData();
  HNDistribute(gf);
}

/*-------------------------------------------------------*/

void StdSolver::AssembleMatrix(const VectorInterface& gu, double d)
{
  _ca.start();
  assert(GetMatrix());

  const GlobalVector& u = GetGV(gu);
  HNAverage(gu);
  HNAverageData();

  GetDiscretization()->Matrix(*GetMatrix(),u,*GetProblemDescriptor()->GetEquation(),d);

  // Face
  if (GetFaceDiscretization())
    {
      const FaceEquation* FEQ = GetProblemDescriptor()->GetFaceEquation();
      assert(FEQ);
      GetFaceDiscretization()->FaceMatrix(*GetMatrix(),u,*FEQ,d);
    }
  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
  {
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetDiscretization()->BoundaryMatrix(*GetMatrix(),u,BM->GetBoundaryEquationColors(),*BE,d);
  }
  PeriodicMatrix();
  DirichletMatrix();
  HNZero(gu);
  HNZeroData();

  _ca.stop();
}

/*-------------------------------------------------------*/

void StdSolver::DirichletMatrix() const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntSet& Colors = BM->GetDirichletDataColors();
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM->GetDirichletDataComponents(col);
      GetDiscretization()->StrongDirichletMatrix(*GetMatrix(), col, comp);
    }
}

/* -------------------------------------------------------*/

void StdSolver::DirichletMatrixOnlyRow() const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntSet& Colors = BM->GetDirichletDataColors();
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM->GetDirichletDataComponents(col);
      GetDiscretization()->StrongDirichletMatrixOnlyRow(*GetMatrix(), col, comp);
    }
}

/* -------------------------------------------------------*/

void StdSolver::PeriodicMatrix() const
{
  /*-------------------------------------------------------
  | Modifiziert die Systemmatrix, um den periodischen
  | Raendern Rechnung zu tragen.
  | Vgl. DirichletMatrix bzw. DirichletMatrixOnlyRow.
  | Ruft dazu Matrix->periodic() auf.
  -------------------------------------------------------*/

  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntVector& iv_PeriodicColors = BM->GetPeriodicDataColors();

  const MeshInterface* p_mesh        = GetMesh();
  const GascoigneMesh* GMP           = dynamic_cast<const GascoigneMesh*>(p_mesh);
  assert(GMP);

  map<int,map<int,int> > mm_PeriodicPairs = GMP->GetBoundaryIndexHandler().GetPeriodicPairs();

  for (IntVector::const_iterator p_col=iv_PeriodicColors.begin(); p_col!=iv_PeriodicColors.end();)
    {
      int col = *p_col++;
      *p_col++;
      
      const IntVector iv_PeriodicComponents = BM->GetPeriodicDataComponents(col);

      GetMatrix()->periodic(mm_PeriodicPairs[col],iv_PeriodicComponents);
    }
}

/* -------------------------------------------------------*/

void StdSolver::ComputeIlu() const
{
#ifdef __WITH_UMFPACK__
  if(_directsolver&&_useUMFPACK)
    {
      _cs.start();
      UmfIlu* UM = dynamic_cast<UmfIlu*>(GetIlu());
      assert(UM);
//       if(PrimalSolve==0) return;
      UM->Factorize();
      _cs.stop();
    }
  else
#endif
    if (GetSolverData().GetLinearSmooth()=="ilu")
    {
      _ci.start();
      IntVector perm(GetIlu()->n());
      iota(perm.begin(),perm.end(),0);
      GetIlu()->ConstructStructure(perm,*GetMatrix());
      GetIlu()->zero();
      GetIlu()->copy_entries(GetMatrix());
      int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
      modify_ilu(*GetIlu(),ncomp);
      GetIlu()->compute_ilu();
      _ci.stop();
    }
}

/* -------------------------------------------------------*/

void StdSolver::ComputeIlu(const VectorInterface& gu) const
{
#ifdef __WITH_UMFPACK__
  if(_directsolver&&_useUMFPACK)
    {
      _cs.start();
      UmfIlu* UM = dynamic_cast<UmfIlu*>(GetIlu());
      assert(UM);
//       if(PrimalSolve==0) return;
      UM->Factorize();
      _cs.stop();
    }
  else
#endif
    if (GetSolverData().GetLinearSmooth()=="ilu")
    {
      int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
      _ci.start();
      PermutateIlu(gu);
      GetIlu()->zero();
      GetIlu()->copy_entries(GetMatrix());
      modify_ilu(*GetIlu(),ncomp);
      GetIlu()->compute_ilu();
      _ci.stop();
    }
}

/*-------------------------------------------------------*/

void StdSolver::modify_ilu(IluInterface& I,int ncomp) const 
{
  if(GetSolverData().GetIluModify().size()==0) return;
  if( GetSolverData().GetIluModify().size()!=ncomp ) {
    cerr << "ERROR: GetSolverData().GetIluModify().size()="<< GetSolverData().GetIluModify().size() << " and ";
    cerr << "ncomp="<< ncomp << endl; 
    abort();
    // assert(GetSolverData().GetIluModify().size()==ncomp);
  }

  for(int c=0;c<ncomp;c++)
    {
      double s = GetSolverData().GetIluModify(c);
      I.modify(c,s);
    }
}

/* -------------------------------------------------------*/

void StdSolver::PermutateIlu(const VectorInterface& gu) const
{
  const GlobalVector& u = GetGV(gu);
  
  #ifdef __WITH_THREADS__
  if(__n_threads >1)
  {
    assert(__n_threads == __thread_domain2node.size());
    
    if (GetSolverData().GetIluSort()=="cuthillmckee")
    {
      std::vector<IntVector> perm(__n_threads);
      for(int d=0; d< __n_threads ; d++)
      {
	int n = __thread_domain2node[d].size();
	perm[d].resize(n); 
	iota(perm[d].begin(),perm[d].end(),0);
	
	CuthillMcKee    cmc(GetMatrix()->GetStencil());
	cmc.Permutate      (perm[d],__thread_domain2node[d],__thread_node2domain,d);   
	
      } 
      //Wie das?
      ThreadIlu* TIlu = dynamic_cast<ThreadIlu*>(GetIlu());
      assert(TIlu);
      TIlu->ConstructStructure(perm,*GetMatrix(),__thread_domain2node);
    }
    else
    {
      std::cerr<<"In StdSolver::PermutateIlu: IluSort "<<GetSolverData().GetIluSort()<<" is not available with threads."<<std::endl;
      abort();
    }
    //end of the thread section return to allow for special case if only one thread is used
    return;
  }
#endif
  //keine Threads oder nur einer zur Verfuegung
  int n = GetMatrix()->GetStencil()->n();
  IntVector perm(n);

  iota(perm.begin(),perm.end(),0);
  if (GetSolverData().GetIluSort()=="cuthillmckee")
    {
      CuthillMcKee    cmc(GetMatrix()->GetStencil());
      cmc.Permutate      (perm);
    }
  else if (GetSolverData().GetIluSort()=="streamdirection")
    {
      assert(GetProblemDescriptor()->GetEquation());
      assert(GetSolverData().GetStreamDirection().size()<=GetProblemDescriptor()->GetEquation()->GetNcomp());
      StreamDirection sd (GetMesh(),GetMatrix()->GetStencil(),u);
      sd.Permutate       (perm,GetSolverData().GetStreamDirection());
    }
  else if (GetSolverData().GetIluSort()=="vectordirection")
    {
      VecDirection vd (GetMesh());
      vd.Permutate    (perm,GetSolverData().GetVectorDirection());
    }
  GetIlu()->ConstructStructure(perm,*GetMatrix());
}

/* -------------------------------------------------------*/

void StdSolver::Visu(const string& name, const VectorInterface& gu, int i) const
{
  if(gu.GetType()=="node")
  {
    PointVisu(name,GetGV(gu),i);
  }
  else if(gu.GetType()=="cell")
  {
    CellVisu(name,GetGV(gu),i);
  }
  else
  {
    cerr << "No such vector type: " << gu.GetType() << endl;
    abort();
  }
}

/* -------------------------------------------------------*/

GascoigneVisualization* StdSolver::NewGascoigneVisualization() const {
  // by virtualizing the instantiation of the GascoigneVisualization object two things can be achieved
  // * setting the time in the overloaded instantiation in StdTimeSolver
  // * instantiating a derived class of GascoigneVisualization = good for transforming the output grid
  GascoigneVisualization* p_gv = new GascoigneVisualization;
  return p_gv;
}

/* -------------------------------------------------------*/

void StdSolver::PointVisu(const string& name, const GlobalVector& u, int i) const
{
  GetDiscretization()->HNAverage(const_cast<GlobalVector&>(u)); 

  GascoigneVisualization* p_gv   = NewGascoigneVisualization();
  GascoigneVisualization& Visu   = *p_gv;

  Visu.SetMesh(GetMesh());  

  const ComponentInformation*  CI = GetProblemDescriptor()->GetComponentInformation();
  if(CI)
  {
    Visu.AddPointVector(CI,&u);
  }
  else
  {
    Visu.AddPointVector(&u);
  }

  Visu.read_parameters(_paramfile);
  Visu.set_name(name);
  Visu.step(i);
  Visu.write();

  delete p_gv;
  GetDiscretization()->HNZero(const_cast<GlobalVector&>(u)); 
}

/* -------------------------------------------------------*/

void StdSolver::CellVisu(const string& name, const GlobalVector& u, int i) const
{
  GascoigneVisualization* p_gv   = NewGascoigneVisualization();
  GascoigneVisualization& Visu   = *p_gv;

  Visu.SetMesh(GetMesh());

  const ComponentInformation*  CI = GetProblemDescriptor()->GetComponentInformation();
  if(CI)
  {
    Visu.AddCellVector(CI,&u);
  }
  else
  {
    Visu.AddCellVector(&u);
  }

  Visu.read_parameters(_paramfile);
  Visu.set_name(name);
  Visu.step(i);
  Visu.write();

  delete p_gv;
}

/* -------------------------------------------------------*/

void StdSolver::VisuGrid(const string& name, int i) const
{
  assert(GetMesh());
  
  if (GetMesh()->dimension()==2)
    {
      VisuEPS eps(_paramfile);
      //  eps.SetOption(VisuEPS::LINEWIDTH,0.1);
      if(_discname[1]=='2')
      {
        eps.SetOption(VisuEPS::WRITE_PATCH,1);
      }
      eps.SetMesh(*GetMesh());
      eps.WriteGrid(name,i);
    }
}

/*-------------------------------------------------------*/

void StdSolver::Read(VectorInterface& gu, const string& filename) const
{
  GlobalVector& u = GetGV(gu);
  u.zero();
  ReadBackUp(u,filename);
}

/*-------------------------------------------------------*/

void StdSolver::Write(const VectorInterface& gu, const string& filename) const
{
  const GlobalVector& u = GetGV(gu);
  WriteBackUp(u,filename);
}

/*-------------------------------------------------------*/

void StdSolver::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
{
  GetDiscretization()->ConstructInterpolator(I,MT);
}

/* -------------------------------------------------------*/

DoubleVector StdSolver::IntegrateSolutionVector(const VectorInterface& gu) const
{
  HNAverage(gu);
  DoubleVector dst = GetPfilter().IntegrateVector(GetGV(gu));
  HNZero(gu);
  return dst;
}

/* -------------------------------------------------------*/

void StdSolver::SubtractMean(VectorInterface& gx) const
{
  GlobalVector& x = GetGV(gx);
  // In each nonlinear step: applied to Newton correction,
  // in each smoothing step
  //
  if (GetPfilter().Active())
    {
      GetDiscretization()->HNZeroCheck(x);
      GetPfilter().SubtractMean(x);
      HNZero(gx);
    }
}

/* -------------------------------------------------------*/

void StdSolver::SubtractMeanAlgebraic(VectorInterface& gx) const
{
  GlobalVector& x = GetGV(gx);
  
  // applies to residuals
  if (GetPfilter().Active())
    {
      GetDiscretization()->HNZeroCheck(x);
      GetPfilter().SubtractMeanAlgebraic(x);
      HNZero(gx);
    }
}

/*---------------------------------------------------*/

void StdSolver::DeleteVector(VectorInterface& p) const
{
  _NGVA.Delete(p);
}

/*-----------------------------------------*/

void StdSolver::Equ(VectorInterface& dst, double s, const VectorInterface& src) const
{
  GetGV(dst).equ(s,GetGV(src));
}

/*-----------------------------------------*/

void StdSolver::Add(VectorInterface& dst, double s, const VectorInterface& src) const
{
  GetGV(dst).add(s,GetGV(src));
}

/*-----------------------------------------*/

void StdSolver::SAdd(double s1,VectorInterface& dst, double s2, const VectorInterface& src) const
{
  GetGV(dst).sadd(s1,s2,GetGV(src));
}

/*-----------------------------------------*/

double StdSolver::Norm(const VectorInterface& dst) const
{
  return GetGV(dst).norm();
}

/*-----------------------------------------*/

double StdSolver::ScalarProduct(const VectorInterface& y, const VectorInterface& x) const
{
  return GetGV(y)*GetGV(x);
}

/*---------------------------------------------------*/

void StdSolver::AssembleDualMatrix(const VectorInterface& gu, double d)
{
  _ca.start();

  MatrixInterface* M = GetMatrix();

  assert(M);

  HNAverage(gu);

  const Equation& EQ = *GetProblemDescriptor()->GetEquation();
  M->zero();
  GetDiscretization()->Matrix(*M,GetGV(gu),EQ,d);
  M->transpose();

  // PeriodicMatrix() hier nicht getestet!
  PeriodicMatrix();
  DirichletMatrixOnlyRow();
  HNZero(gu);

  _ca.stop();
}

/*---------------------------------------------------*/

void StdSolver::RhsCurve(VectorInterface &f, const Curve &C,int comp,int N) const
{
  HNAverageData();

  GetDiscretization()->RhsCurve(GetGV(f),C,comp,N);

  HNZeroData();
  HNDistribute(f);
}
/*--------------------------------------------------------*/

double StdSolver::ScalarProductWithFluctuations(DoubleVector& eta, const VectorInterface& gf, 
						const VectorInterface& gz) const
{
  const GlobalVector& f = GetGV(gf);
  const GlobalVector& z = GetGV(gz);

  GlobalVector dz(f.ncomp(),f.n());

  dz.zero();
  Pi pi;
  pi.Init(GetMesh());
  pi.vmult(dz,z);

  for(int i=0; i<z.n(); i++)
    {
      for (int c=0; c<z.ncomp(); c++)
	{
	  eta[i] += fabs(f(i,c)*dz(i,c));
	}
    } 
  return dz * f;
}

}
