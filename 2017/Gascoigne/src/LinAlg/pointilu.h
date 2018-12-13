/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#ifndef  __PointIlu_h
#define  __PointIlu_h

#include  "iluinterface.h"
#include  "simpleilu.h"
#include  "sparsestructureadaptor.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments PointIlu

///
///
/////////////////////////////////////////////

class PointIlu : public SimpleIlu, virtual public IluInterface
{
protected:

  int _ncomp;
  SparseStructureAdaptor* SSAP;

public:

//
///  Constructor 
//
    
    PointIlu(int ncomp, std::string type);
    ~PointIlu();

    std::string GetName() const {return "PointIlu";}
    
    void ReInit(const SparseStructureInterface* S);
    
    int   n()          const { return GetStencil()->n();};
    void zero() {
      SimpleIlu::zero();
    }
    
    void ConstructStructure(const IntVector& perm, const MatrixInterface& A);
    void modify(int c, double s);
    void copy_entries(const MatrixInterface*  A) {
      SimpleIlu::copy_entries(A);
    }
    void compute_ilu() {
      SimpleIlu::compute_ilu();
    }
    void solve(GlobalVector& x) const {
      SimpleIlu::solve(x);
    }
    void solve_transpose(GlobalVector& x) const {
      SimpleIlu::solve_transpose(x);
    }
};
}

#endif
