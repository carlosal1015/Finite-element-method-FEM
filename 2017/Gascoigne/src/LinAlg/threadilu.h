/**
*
* Copyright (C) 2009, 2011 by the Gascoigne 3D authors
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


/*----------------------------   threadilu.h     ---------------------------*/
/*      $Id$                 */
#ifndef __threadilu_H
#define __threadilu_H
/*----------------------------   threadilu.h     ---------------------------*/


#include  "threadsparseblockilu.h"
#include  <vector>

/*-------------------------------------------------------------*/
#ifdef __WITH_THREADS__
namespace Gascoigne
{
  class ThreadIlu : public IluInterface
  {
    private:
    
    protected:

    std::vector<ThreadSparseBlockIluInterface*>     __local_ilu;

    int                            __n_threads, __n, __ncomp;
    std::vector<std::vector<int> > __permutations;
    std::vector<std::vector<int> > __domains2nodes;
    nvector<double>                __weights; 

    mutable std::vector<GlobalVector>      __local_vectors;

  public:
    ThreadIlu(int ncomp);
    virtual ~ThreadIlu();

    virtual int   n() const             {
      std::cerr << "\"ThreadILU::n\" not written!" << std::endl;
      abort();
//       return __n;
    }
    virtual std::string GetName() const {return "Thread-Ilu";}
    
    virtual void ReInit(const SparseStructureInterface* A);
    
    virtual void ConstructStructure(const  std::vector<IntVector>& perm, const MatrixInterface& A, const std::vector<std::vector<int> >& domains2nodes);

    virtual void ConstructStructure(const Gascoigne::IntVector&, const Gascoigne::MatrixInterface&)
    {
      std::cerr <<"\"ThreadIlu::ConstructStructure(const Gascoigne::IntVector&, const Gascoigne::MatrixInterface&)\" has no meaning, please use !" << std::endl;
      std::cerr <<"\"ThreadIlu::ConstructStructure(const  std::vector<std::vector<int> >& perm, const MatrixInterface& A, const std::vector<std::vector<int> > domains2nodes)\" instead!" << std::endl;
     abort();
    }

    virtual void modify(int c, double s);
    virtual void zero() ;
    virtual void compute_ilu ();
    virtual void copy_entries(const MatrixInterface* A);
    
    virtual void solve(GlobalVector& x) const;
  };
}
//endof __WITH_THREADS__
#endif

/*----------------------------   threadilu.h     ---------------------------*/
/* end of #ifndef __threadilu_H */
#endif
/*----------------------------   threadilu.h     ---------------------------*/
