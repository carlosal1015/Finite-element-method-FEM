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


#ifndef  __FiniteElementWithSecond_h
#define  __FiniteElementWithSecond_h

#include  "finiteelement.h"

namespace Gascoigne
{
template<int DIM, int BDIM, class TRAFO, class BASE>
class FiniteElementWithSecond : public FiniteElement<DIM,BDIM,TRAFO,BASE>
{
  protected:
    
    typedef  FemInterface::Matrix   Matrix;
    
    mutable nvector<Matrix> hesse;
    
  public:

    void ComputeHesse(const Vertex2d& xi) const;
    void ComputeHesse(const Vertex3d& xi) const;

    std::string GetName() const {return "FiniteElementWithSecond";}
    
    FiniteElementWithSecond();
    
    void point(const Vertex<DIM>& v) const
    {
      FiniteElement<DIM,BDIM,TRAFO,BASE>::point(v);
      ComputeHesse(v);
    }

    void  init_test_functions(TestFunction& Phi, double w, int i) const
    {
      FiniteElement<DIM,BDIM,TRAFO,BASE>::init_test_functions(Phi,w,i);
      init_test_hesse(Phi,w,i);
    }
    void init_test_hesse(TestFunction& N, double w, int i) const;
};
}

#endif

