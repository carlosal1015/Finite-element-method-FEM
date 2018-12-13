/**
*
* Copyright (C) 2004, 2005, 2011 by the Gascoigne 3D authors
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


#ifndef __sparseblockilu_h
#define __sparseblockilu_h

#include  "sparseblockmatrix.h"
#include  "iluinterface.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
template<class B>
class SparseBlockIlu: public virtual IluInterface, public SparseBlockMatrix<B>
{
protected:

  nvector<int>          p,q;
  GlobalVector*   yp;

  void backward() const;
  void forward () const;
  virtual void hin(const GlobalVector& x) const;
  virtual void her(GlobalVector& x) const;

  int   n()          const { return SparseBlockMatrix<B>::US.n();};
  const int&  start(int i) const { return SparseBlockMatrix<B>::US.start(i); }; 
  const int&  stop (int i) const { return SparseBlockMatrix<B>::US.stop (i); }; 
  const int&  col(int pos) const { return SparseBlockMatrix<B>::US.col(pos); };
  const int&  diag(int i)  const { return SparseBlockMatrix<B>::US.diag(i); }; 

  public:

  SparseBlockIlu<B>();
  SparseBlockIlu<B>(const SparseBlockIlu<B>& I);
  ~SparseBlockIlu();

  string GetName() const {return "SparseBlockIlu";}
  
  nvector<int>&       GetP() {return p;}
  nvector<int>&       GetQ() {return q;}
  const nvector<int>& GetP() const {return p;}
  const nvector<int>& GetQ() const {return q;}

  void modify(int c, double s);
  void zero() { SparseBlockMatrix<B>::zero(); }

  void compute_ilu ();
  void ReInit      (const SparseStructureInterface* SI);
  void ConstructStructure(const nvector<int>& perm, const MatrixInterface& A);
  void copy_entries(const MatrixInterface* A);
  void solve       (GlobalVector& x) const;
  void solvetrans  (GlobalVector& x) const {
    std::cerr << "\"SparseBlockIlu::solvetrans\" not written!" << std::endl;
    abort();
  }
  ostream& Write(ostream &s) const;
};
}

#endif
