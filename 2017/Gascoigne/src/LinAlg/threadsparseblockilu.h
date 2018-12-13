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


#ifndef __threadsparseblockilu_h
#define __threadsparseblockilu_h

#include  "sparseblockilu.h"


#ifdef __OLDCOMPILER__
  #include  <hash_map>
  #define HASHMAP  hash_map
#else
#ifdef __NEWER_THAN_GCC_4_2__
  #include <tr1/unordered_map>
  #define HASHMAP   std::tr1::unordered_map
#else
  #include  <ext/hash_map>
  #define HASHMAP  __gnu_cxx::hash_map
#endif
#endif


/*-------------------------------------------------------------*/

namespace Gascoigne
{

  class ThreadSparseBlockIluInterface : public IluInterface
  {
  public:
    virtual void ConstructStructure(const nvector<int>& perm, const MatrixInterface& A,
				    const nvector<int>& nodes_of_domain)
    {
      std::cerr << "ThreadSparseBlockIluInterface::ConstructStructure" << endl;
      abort();
    }
    
    virtual const nvector<int>&     GetP() const=0;
    virtual const HASHMAP<int,int>& GetQ() const=0;
//    virtual int   n()                      const=0;
  };

  // --------------------------------------------------

template<class B>
  class ThreadSparseBlockIlu  :  public SparseBlockMatrix<B>, public ThreadSparseBlockIluInterface
  {
protected:
    
    nvector<int>     p;
    HASHMAP<int,int> q;

    void backward(GlobalVector& x) const;
    void forward (GlobalVector& x) const;

    int   n()          const { return SparseBlockMatrix<B>::US.n();};
    const int&  start(int i) const { return SparseBlockMatrix<B>::US.start(i); }; 
    const int&  stop (int i) const { return SparseBlockMatrix<B>::US.stop (i); }; 
    const int&  col(int pos) const { return SparseBlockMatrix<B>::US.col(pos); };
    const int&  diag(int i)  const { return SparseBlockMatrix<B>::US.diag(i); }; 
    
  public:

  ThreadSparseBlockIlu<B>();
  ThreadSparseBlockIlu<B>(const ThreadSparseBlockIlu<B>& I);
  ~ThreadSparseBlockIlu();

  string GetName() const {return "ThreadSparseBlockIlu";}
  
  const nvector<int>&     GetP() const {return p;}
  const HASHMAP<int,int>& GetQ() const {return q;}

  void modify(int c, double s);
  void zero() { SparseBlockMatrix<B>::zero(); }

  void compute_ilu ();
  void ReInit      (const SparseStructureInterface* SI);
  void ConstructStructure(const nvector<int>& perm, const MatrixInterface& A)
  {abort();}
  void ConstructStructure(const nvector<int>& perm, const MatrixInterface& A,
			  const nvector<int>& nodes_of_domain);

  void copy_entries(const MatrixInterface* A);
  void solve       (GlobalVector& x) const;
  void solvetrans  (GlobalVector& x) const {
    std::cerr << "\"ThreadSparseIlu::solvetrans\" not written!" << std::endl;
    abort();
  }
  ostream& Write(ostream &s) const;
};
}

#undef HASHMAP

#endif
