/**
*
* Copyright (C) 2004, 2005, 2010 by the Gascoigne 3D authors
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


#ifndef __FiniteElement_h
#define __FiniteElement_h

#include  "feminterface.h"

/*-----------------------------------------------------*/

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  FE based on Transformation (TRAFO) and Referenzelement (BASE)

///
///
/////////////////////////////////////////////

template<int DIM, int BDIM, class TRAFO, class BASE>
class FiniteElement : public FemInterface
{
 protected:

  TRAFO                   T;
  BASE                    B;
  mutable std::vector<Vertex<DIM> >    grad;
  mutable double                  det;

  virtual void ComputeGrad() const;

 public:

  FiniteElement();

  std::string GetName() const {return "FiniteElement";}

  int    n()          const { return B.n(); }
  double N   (int i)  const { return B.phi(i); }
  double N_x (int i)  const { return grad[i].x(); }
  double N_y (int i)  const { return grad[i].y(); }
  double N_z (int i)  const { return grad[i].z(); }
  double J()          const { return det;}
  double G()          const { return T.G();}

  void x     (Vertex<DIM>& v) const { v = T.x();}
  void normal(Vertex<DIM>& v) const { v = T.normal();};

  void point(const Vertex<DIM>&) const;
  void point_boundary(int ie, const Vertex<BDIM>& s1) const;
  /// depreciated
  void ReInit(const Matrix& M) const { T.ReInit(M); }

  void  init_test_functions(TestFunction& Phi, double w, int i) const;

  void Anisotropy(DoubleMatrix& A) const;

  void GetCoordinates(DoubleMatrix& A) const { T.GetCoordinates(A);}
};
}

/*-----------------------------------------------------*/

#endif
