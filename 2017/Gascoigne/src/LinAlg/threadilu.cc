/**
*
* Copyright (C) 2009 by the Gascoigne 3D authors
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


#include "threadilu.h"
#include "sparsestructure.h"
#include "threadsparseblockilu.h"
#include "fmatrixblock.h"
#ifdef __WITH_THREADS__
#include <omp.h>

using namespace std;

namespace Gascoigne
{

  ThreadIlu::ThreadIlu(int ncomp) : __ncomp(ncomp)
  {}

  ThreadIlu::~ThreadIlu()
  {
    for (int i=0;i<__local_ilu.size();++i)
      if (__local_ilu[i]) delete __local_ilu[i];
    __local_ilu.clear();
  }


  // --------------------------------------------------

  void ThreadIlu::ReInit(const SparseStructureInterface* A)
  {
  }

   // --------------------------------------------------

  void ThreadIlu::ConstructStructure(const std::vector<IntVector>& perm, const MatrixInterface& MAP, const std::vector<std::vector<int> >& domains2nodes)
  {

    __n_threads = perm.size();
    
    __permutations.resize(__n_threads);
#pragma omp parallel for schedule(static,1)
    for (int i=0;i<__n_threads;++i)
      __permutations[i] = perm[i];

    __domains2nodes.resize(__n_threads);
#pragma omp parallel for schedule(static,1)
    for (int i=0;i<__n_threads;++i)
      __domains2nodes[i] = domains2nodes[i];

    __local_vectors.resize(__n_threads);
#pragma omp parallel for schedule(static,1)
    for (int i=0;i<__n_threads;++i)
    {
      __local_vectors[i].ncomp() = __ncomp;
      __local_vectors[i].resize(perm[i].size());
    }

    
    __weights.resize(MAP.GetStencil()->n());
    __weights.zero();
    
    for (int p=0;p<perm.size();++p)
    {
#pragma omp parallel for 
      for (int i=0;i<perm[p].size();++i)
	__weights[domains2nodes[p][perm[p][i]]] += 1.0;
    }

#pragma omp parallel for
   for (int i=0;i<__weights.size();++i)
    {
      assert(__weights[i]>0);
      __weights[i] = 1.0 / __weights[i];
    }
    

    // new ilu
    for (int i=0;i<__local_ilu.size();++i)
    {
      assert(__local_ilu[i]);
      delete __local_ilu[i];
    }
    __local_ilu.resize(__n_threads);
#pragma omp parallel for schedule(static,1)
   for (int i=0;i<__n_threads;++i)
    {
      if (__ncomp == 1) __local_ilu[i] = new ThreadSparseBlockIlu<FMatrixBlock<1> >;
      if (__ncomp == 2) __local_ilu[i] = new ThreadSparseBlockIlu<FMatrixBlock<2> >;
      if (__ncomp == 3) __local_ilu[i] = new ThreadSparseBlockIlu<FMatrixBlock<3> >;
      if (__ncomp == 4) __local_ilu[i] = new ThreadSparseBlockIlu<FMatrixBlock<4> >;
      if (__ncomp == 5) __local_ilu[i] = new ThreadSparseBlockIlu<FMatrixBlock<5> >;
      if (__ncomp == 6) __local_ilu[i] = new ThreadSparseBlockIlu<FMatrixBlock<6> >;
      if (__ncomp == 7) __local_ilu[i] = new ThreadSparseBlockIlu<FMatrixBlock<7> >;
      if (__ncomp == 8) __local_ilu[i] = new ThreadSparseBlockIlu<FMatrixBlock<8> >;
    }
      
    
    // construct stencils
#pragma omp parallel for schedule(static,1)
   for (int p=0;p<__n_threads;++p)
    {
      __local_ilu[p]->ConstructStructure(perm[p],MAP,domains2nodes[p]);
    }
  }


  // --------------------------------------------------

  void ThreadIlu::zero()
  {
#pragma omp parallel for schedule(static,1)
    for (int p=0;p<__n_threads;++p)
      __local_ilu[p]->zero();
  }
  
  // --------------------------------------------------
  
  void ThreadIlu::modify(int c, double s)
  {
#pragma omp parallel for schedule(static,1)
    for (int p=0;p<__n_threads;++p)
      __local_ilu[p]->modify(c,s);
  }
  
  // --------------------------------------------------
  
  void ThreadIlu::copy_entries(const MatrixInterface* A) 
  {
#pragma omp parallel for schedule(static,1)
    for (int p=0;p<__n_threads;++p)
      __local_ilu[p]->copy_entries(A);
  }

  // --------------------------------------------------
  
  void ThreadIlu::compute_ilu()
  {
#pragma omp parallel for schedule(static,1)
   for (int p=0;p<__n_threads;++p)
      __local_ilu[p]->compute_ilu();
  }

  // --------------------------------------------------

  void ThreadIlu::solve(GlobalVector& x) const 
  {
 #pragma omp parallel for schedule(static,1)
   for (int p=0;p<__n_threads;++p)
    {
      assert(__local_vectors[p].n() == __local_ilu[p]->n());
      for (int i=0;i<__local_ilu[p]->n();++i)
	__local_vectors[p].equ_node(i, 1.0, __local_ilu[p]->GetP()[i], x);
      __local_ilu[p]->solve(__local_vectors[p]);
    }

    x.zero();
    for (int p=0;p<__n_threads;++p)
    {
      // hier ueberschneiden sich die vektoren am interface
#pragma omp parallel for
     for (int i=0;i<__local_ilu[p]->n();++i)
      {
	int i_global = __local_ilu[p]->GetP()[i];
	double weight = __weights[i_global];
	x.add_node(i_global, weight , i , __local_vectors[p]);
      }
    }
    
    
  }
 


}
//endof __WITH_THREADS__
#endif
