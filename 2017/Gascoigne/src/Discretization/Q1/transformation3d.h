/**
*
* Copyright (C) 2004, 2005, 2006, 2010 by the Gascoigne 3D authors
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


#ifndef __Transformation3d_h
#define __Transformation3d_h

#include  "nmatrix.h"
#include  "vertex.h"

/*-----------------------------------------------------*/

namespace Gascoigne
{
template<class BASE>
class Transformation3d
{
 protected:

  typedef nmatrix<double>    Matrix;

  BASE            B;
  mutable Matrix  X;
  mutable Matrix  dt, dti;

  inline void ComputeDT() const;
  	    // second derivatives tensor
  const nvector<Matrix>  ComputeDDT (const Vertex3d& xi) const;

 public:

  Transformation3d();

  const Matrix& DT () const {return dt ;}
  const Matrix& DTI() const {return dti;}

  	    // inverse of second derivatives tensor
  const nvector<Matrix>  DDTI(const Vertex3d& xi) const;
  
  inline double        J     () const;
  inline double        G     () const;
  inline Vertex3d      x     () const;
  inline Vertex3d      normal() const;
  inline void  init          (const Matrix& M) {X=M;}
  inline void  ReInit          (const Matrix& M) const {X=M;}
  inline void  point         (const Vertex3d& xi) const;
  inline void  point_boundary(int ie, const Vertex2d& s) const;
  inline void  GetCoordinates(Matrix& A) const { A.equ(1.,X);}
};

/*-----------------------------------------------------*/

template<class BASE>
inline Transformation3d<BASE>::Transformation3d() : B()
{
  X  .memory(3,B.n());
  dt .memory(3,3);
  dti.memory(3,3);
}

/*-----------------------------------------------------*/

template<class BASE>
inline Vertex3d  Transformation3d<BASE>::x() const 
{
  Vertex3d xp;
  for(int i=0;i<B.n();i++)
    {
      xp.x() += X(0,i) * B.phi(i);
      xp.y() += X(1,i) * B.phi(i);
      xp.z() += X(2,i) * B.phi(i);
    }
  return xp;
}

/*-----------------------------------------------------*/

template<class BASE>
inline Vertex3d  Transformation3d<BASE>::normal() const 
{
  Vertex3d xn;
  dti.mult(xn,*B.normal3d());
  double xx = sqrt(xn*xn);
  xn /= xx;
  return xn;
}

/*-----------------------------------------------------*/

template<class BASE>
inline void  Transformation3d<BASE>::ComputeDT() const
{
  dt.zero();
  for(int i=0;i<B.n();i++)
    {
      dt(0,0) += X(0,i) * B.phi_x(i);
      dt(0,1) += X(0,i) * B.phi_y(i);
      dt(0,2) += X(0,i) * B.phi_z(i);
      dt(1,0) += X(1,i) * B.phi_x(i);
      dt(1,1) += X(1,i) * B.phi_y(i);
      dt(1,2) += X(1,i) * B.phi_z(i);
      dt(2,0) += X(2,i) * B.phi_x(i);
      dt(2,1) += X(2,i) * B.phi_y(i);
      dt(2,2) += X(2,i) * B.phi_z(i);
    }
  dti(0,0) = dt(0,0);
  dti(0,1) = dt(1,0);
  dti(0,2) = dt(2,0);
  dti(1,0) = dt(0,1);
  dti(1,1) = dt(1,1);
  dti(1,2) = dt(2,1);
  dti(2,0) = dt(0,2);
  dti(2,1) = dt(1,2);
  dti(2,2) = dt(2,2);

  dti.gauss_jordan();
}


/*-----------------------------------------------------*/

template<class BASE>
inline const nvector<nmatrix<double> > Transformation3d<BASE>::DDTI(const Vertex3d& xi) const
{
  const nvector<nmatrix<double> >& ddt = ComputeDDT(xi);
  
  nvector<nmatrix<double> > ddti(3,nmatrix<double> (3,3));
  Matrix dti_ = dti;
  

  dti_.transpose();
  for (int i=0;i<3;++i)
    {
      nmatrix<double> tmp(3,3);

      ddti[i].zero();
      for (int d=0;d<3;++d)
	ddti[i].add(-dti_(d,i),ddt[d]);

      ddti[i].mmult(tmp,dti_);
      dti_.mmult(ddti[i],tmp);
    }
  return ddti;  
}

/*-----------------------------------------------------*/
template<class BASE>
inline const nvector<nmatrix<double> > Transformation3d<BASE>::ComputeDDT(const Vertex3d& xi) const
{
  const_cast<BASE*> (&B)->point(xi);
  
  nvector<nmatrix<double> > ddt(3,nmatrix<double> (3,3));
  for (int i=0;i<3;++i) ddt[i].zero();

  for (int i=0;i<B.n();++i)
    {
      for (int j=0;j<3;++j)
	{
	  ddt[0](j,0) += X(j,i) * B.phi_xx(i);
	  ddt[0](j,1) += X(j,i) * B.phi_xy(i);
	  ddt[0](j,2) += X(j,i) * B.phi_xz(i);

	  ddt[1](j,0) += X(j,i) * B.phi_xy(i);
	  ddt[1](j,1) += X(j,i) * B.phi_yy(i);
	  ddt[1](j,2) += X(j,i) * B.phi_yz(i);
	  
	  ddt[2](j,0) += X(j,i) * B.phi_xz(i);
	  ddt[2](j,1) += X(j,i) * B.phi_yz(i);
	  ddt[2](j,2) += X(j,i) * B.phi_zz(i);
	}
    }
  return ddt;  
}

/*-----------------------------------------------------*/

template<class BASE>
inline void  Transformation3d<BASE>::point(const Vertex3d& xi) const
{
  B.point(xi);
  ComputeDT();
}

/*-----------------------------------------------------*/

template<class BASE>
inline void  Transformation3d<BASE>::point_boundary(int ie, const Vertex2d& s) const
{
  B.point_boundary(ie,s);
  ComputeDT();
}

/*-----------------------------------------------------*/

template<class BASE>
inline double  Transformation3d<BASE>::J() const  
{
  double h = dt.det();
  if (h<0)
    {
      std::cout << "Transformation3d<BASE>::J()  h = " << h << std::endl;
    }
  return h;
}

/*-----------------------------------------------------*/

template<class BASE>
inline double  Transformation3d<BASE>::G() const  
{
  double d1phi=0,d2phi=0,d12phi=0;
  const fixarray<2,int>& fc = *B.faces();
  for (int i=0;i<3;++i)
    {
      d1phi+=dt(i,fc[0])*dt(i,fc[0]);
      d2phi+=dt(i,fc[1])*dt(i,fc[1]);
      d12phi+=dt(i,fc[0])*dt(i,fc[1]);
    }
  double h = d1phi*d2phi-d12phi*d12phi;
  assert(h>=0);
  return sqrt(h);
}
}

#endif
