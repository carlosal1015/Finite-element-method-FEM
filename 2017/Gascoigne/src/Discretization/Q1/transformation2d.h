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


#ifndef __Transformation2d_h
#define __Transformation2d_h

#include  "nmatrix.h"
#include  "vertex.h"

/*-----------------------------------------------------*/

namespace Gascoigne
{
template<class BASE>
class Transformation2d
{
 protected:

  typedef nmatrix<double>    Matrix;

  BASE            B;
  mutable Matrix  X;
  mutable Matrix  dt, dti;

  inline void ComputeDT() const;

  	    // second derivatives tensor
  const nvector<Matrix>  ComputeDDT (const Vertex2d& xi) const;

 public:

  Transformation2d();

  const Matrix& DT () const {return dt ;}
  const Matrix& DTI() const {return dti;}

	    // inverse of second derivatives tensor
  const nvector<Matrix>  DDTI(const Vertex2d& xi) const;
  
  inline double        J     () const;
  inline double        G     () const;
  inline Vertex2d      x     () const;
  inline Vertex2d      normal() const;
  inline void  init          (const Matrix& M) {X=M;}
  inline void  ReInit          (const Matrix& M) const {X=M;}
  inline void  point         (const Vertex2d& xi) const;
  inline void  point_boundary(int ie, const Vertex1d& s) const;
  inline void  GetCoordinates(Matrix& A) const { A.equ(1.,X);}
};

/*-----------------------------------------------------*/

template<class BASE>
inline Transformation2d<BASE>::Transformation2d() : B()
{
  int n = B.n();
  X.memory(2,n);
  dt .memory(2,2);
  dti.memory(2,2);
}

/*-----------------------------------------------------*/

template<class BASE>
inline Vertex2d  Transformation2d<BASE>::x() const 
{
  Vertex2d xp;
  for(int i=0;i<B.n();i++)
    {
      xp.x() += X(0,i) * B.phi(i);
      xp.y() += X(1,i) * B.phi(i);
    }
  return xp;
}

/*-----------------------------------------------------*/

template<class BASE>
inline Vertex2d  Transformation2d<BASE>::normal() const 
{
  Vertex2d xn;
  dti.mult(xn,*B.normal2d());
  double xx = sqrt(xn*xn);
  xn /= xx;
  return xn;
}

/*-----------------------------------------------------*/

template<class BASE>
inline void  Transformation2d<BASE>::ComputeDT() const
{
  dt.zero();
  for(int i=0;i<B.n();i++)
    {
      dt(0,0) += X(0,i) * B.phi_x(i);
      dt(0,1) += X(0,i) * B.phi_y(i);
      dt(1,0) += X(1,i) * B.phi_x(i);
      dt(1,1) += X(1,i) * B.phi_y(i);
    }
  dti(0,0) = dt(0,0);      
  dti(0,1) = dt(1,0);      
  dti(1,0) = dt(0,1);      
  dti(1,1) = dt(1,1);      
  dti.gauss_jordan();
}

/*-----------------------------------------------------*/

template<class BASE>
inline const nvector<nmatrix<double> > Transformation2d<BASE>::ComputeDDT(const Vertex2d& xi) const
{
  const_cast<BASE*> (&B)->point(xi);
  
  nvector<nmatrix<double> > ddt(2,nmatrix<double> (2,2));
  for (int i=0;i<2;++i) ddt[i].zero();

  for (int i=0;i<B.n();++i)
    {
      ddt[0](0,0) += X(0,i) * B.phi_xx(i);
      ddt[0](0,1) += X(0,i) * B.phi_xy(i);
      ddt[0](1,0) += X(1,i) * B.phi_xx(i);
      ddt[0](1,1) += X(1,i) * B.phi_xy(i);

      ddt[1](0,0) += X(0,i) * B.phi_xy(i);
      ddt[1](0,1) += X(0,i) * B.phi_yy(i);
      ddt[1](1,0) += X(1,i) * B.phi_xy(i);
      ddt[1](1,1) += X(1,i) * B.phi_yy(i);
    }
  return ddt;  
}

/*-----------------------------------------------------*/

template<class BASE>
inline const nvector<nmatrix<double> > Transformation2d<BASE>::DDTI(const Vertex2d& xi) const
{
  const nvector<nmatrix<double> >& ddt = ComputeDDT(xi);
  
  nvector<nmatrix<double> > ddti(2,nmatrix<double> (2,2));
  Matrix dti_ = dti;
  

  dti_.transpose();
  for (int i=0;i<2;++i)
    {
      nmatrix<double> tmp(2,2);

      ddti[i].zero();
      for (int d=0;d<2;++d)
	ddti[i].add(-dti_(d,i),ddt[d]);
      
      ddti[i].mmult(tmp,dti_);
      dti_.mmult(ddti[i],tmp);
    }
  return ddti;  
}

/*-----------------------------------------------------*/

template<class BASE>
inline void  Transformation2d<BASE>::point(const Vertex2d& xi) const
{
  B.point(xi);
  ComputeDT();
}

/*-----------------------------------------------------*/

template<class BASE>
inline void  Transformation2d<BASE>::point_boundary(int ie, const Vertex1d& s) const
{
  B.point_boundary(ie,s);
  ComputeDT();
}

/*-----------------------------------------------------*/

template<class BASE>
inline double  Transformation2d<BASE>::J() const  
{
  return dt(0,0)*dt(1,1)-dt(1,0)*dt(0,1);
}

/*-----------------------------------------------------*/

template<class BASE>
inline double  Transformation2d<BASE>::G() const  
{
  Vertex2d xt;
  dt.mult(xt,*B.tangent2d());
  return xt.norm_l2();
}
}

#endif
