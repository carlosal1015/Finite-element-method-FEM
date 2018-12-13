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


#ifndef  __SimpleIlu_h
#define  __SimpleIlu_h

#include  "simplematrix.h"
#include  "matrixinterface.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments SimpleIlu

///
///
/////////////////////////////////////////////

class SimpleIlu : public SimpleMatrix
{
protected:

  IntVector              p,q;
  mutable DoubleVector   yp;

  void hin(const DoubleVector& y) const;
  void her(DoubleVector& y) const;
  void backward() const;
  void forward () const;
  void backward_transpose() const;
  void forward_transpose () const;

public:

//
///  Constructor 
//

    SimpleIlu() : SimpleMatrix() {}
    
      int   n()          const { return GetStencil()->n();};
    void zero() {SimpleMatrix::zero();}
    void ReInit(int n, int nentries);
    void copy_entries(const MatrixInterface*  A);
    void compute_ilu();
    void solve(DoubleVector& x) const;
    void solve_transpose(DoubleVector& x) const;
};
}

#endif
