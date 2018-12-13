/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#ifndef  __UmfIlu_h
#define  __UmfIlu_h

#ifdef __WITH_UMFPACK__

#include  "iluinterface.h"
#include  "simplematrix.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments UmfIlu
///
///
/////////////////////////////////////////////

class UmfIlu : virtual public IluInterface, public SimpleMatrix
{
private:

  const SimpleMatrix* AP;

protected:

  // fuer umfpack
  double *Control;
  double *Info;
  void *Symbolic, *Numeric ;

public:

  //
  ///  Constructor 
    //
    UmfIlu(const MatrixInterface* A);
    ~UmfIlu();
    
    std::string GetName() const { return "UmfIlu"; }
    
    int   n()          const { return GetStencil()->n();};
    void ReInit(const SparseStructureInterface* SS);

    void copy_entries(const MatrixInterface&  A);
    void ConstructStructure(const IntVector& perm, const MatrixInterface& A);
    void Factorize();
    void Solve(DoubleVector& x, const DoubleVector& b);
    void SolveTranspose(DoubleVector& x, const DoubleVector& b);
};
}

#endif
#endif
