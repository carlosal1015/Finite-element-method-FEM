/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#include  "exactsolution.h"

/*-----------------------------------------*/

namespace Gascoigne
{
double ExactSolution::x(int c, const Vertex2d& v)const
{
  Vertex2d vl(v), vr(v);
  vr.x() += eps; vl.x() -= eps;
  return ((*this)(c,vr)-(*this)(c,vl))/(2.*eps);
}
double ExactSolution::xx(int c, const Vertex2d& v)const
{
  Vertex2d vl(v), vr(v);
  vr.x() += eps; vl.x() -= eps;
  return (-2.*(*this)(c,v)+ (*this)(c,vr)+(*this)(c,vl))/(eps*eps);

  return (x(c,vr)-x(c,vl))/(2.*eps);
}
double ExactSolution::y(int c, const Vertex2d& v)const
{
  Vertex2d vl(v), vr(v);
  vr.y() += eps; vl.y() -= eps;
  return ((*this)(c,vr)-(*this)(c,vl))/(2.*eps);
}
double ExactSolution::yy(int c, const Vertex2d& v)const
{
  Vertex2d vl(v), vr(v);
  vr.y() += eps; vl.y() -= eps;
  return (-2.*(*this)(c,v)+ (*this)(c,vr)+(*this)(c,vl))/(eps*eps);

  return (y(c,vr)-y(c,vl))/(2.*eps);
}
double ExactSolution::yx(int c, const Vertex2d& v)const
{
  Vertex2d vl(v), vr(v);
  vr.y() += eps; vl.y() -= eps;
  return (x(c,vr)-x(c,vl))/(2.*eps);
}
double ExactSolution::xy(int c, const Vertex2d& v)const
{
  return yx(c,v);
}

/*-----------------------------------------------------------------*/

double ExactSolution::x(int c, const Vertex3d& v)const
{
  Vertex3d vl(v), vr(v);
  vr.x() += eps; vl.x() -= eps;
  return ((*this)(c,vr)-(*this)(c,vl))/(2.*eps);
}

double ExactSolution::y(int c, const Vertex3d& v)const
{
  Vertex3d vl(v), vr(v);
  vr.y() += eps; vl.y() -= eps;
  return ((*this)(c,vr)-(*this)(c,vl))/(2.*eps);
}

double ExactSolution::z(int c, const Vertex3d& v)const
{
  Vertex3d vl(v), vr(v);
  vr.z() += eps; vl.z() -= eps;
  return ((*this)(c,vr)-(*this)(c,vl))/(2.*eps);
}

double ExactSolution::xx(int c, const Vertex3d& v)const
{
  Vertex3d vl(v), vr(v);
  vr.x() += eps; vl.x() -= eps;
  return (-2.*(*this)(c,v)+ (*this)(c,vr)+(*this)(c,vl))/(eps*eps);
}

double ExactSolution::yy(int c, const Vertex3d& v)const
{
  Vertex3d vl(v), vr(v);
  vr.y() += eps; vl.y() -= eps;
  return (-2.*(*this)(c,v)+ (*this)(c,vr)+(*this)(c,vl))/(eps*eps);
}

double ExactSolution::zz(int c, const Vertex3d& v)const
{
  Vertex3d vl(v), vr(v);
  vr.z() += eps; vl.z() -= eps;
  return (-2.*(*this)(c,v)+ (*this)(c,vr)+(*this)(c,vl))/(eps*eps);
}

double ExactSolution::xy(int c, const Vertex3d& v)const
{
  Vertex3d vl(v), vr(v);
  vr.y() += eps; vl.y() -= eps;
  return (x(c,vr)-x(c,vl))/(2.*eps);
}

double ExactSolution::xz(int c, const Vertex3d& v)const
{
  Vertex3d vl(v), vr(v);
  vr.z() += eps; vl.z() -= eps;
  return (x(c,vr)-x(c,vl))/(2.*eps);
}

double ExactSolution::yz(int c, const Vertex3d& v)const
{
  Vertex3d vl(v), vr(v);
  vr.z() += eps; vl.z() -= eps;
  return (y(c,vr)-y(c,vl))/(2.*eps);
}
}
