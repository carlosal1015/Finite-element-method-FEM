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


#ifndef  __ConstantRightHandSide_h
#define  __ConstantRightHandSide_h

#include  "domainrighthandside.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class OneRightHandSideData : public DomainRightHandSide
{
  protected :

    int ncomp;

  public  : 

    OneRightHandSideData(int n) : DomainRightHandSide(), ncomp(n) {}
    std::string GetName() const { return "one" ;} 
    int GetNcomp() const {return ncomp;}
    double operator()(int c, const Vertex2d& v)const {return 1.;}
    double operator()(int c, const Vertex3d& v)const {return 1.;}
};

/*-----------------------------------------*/

class ConstantRightHandSideData : public DomainRightHandSide
{
protected :
  int     _comp,_ncomp;
  double  _d;

public  : 
  ConstantRightHandSideData(const std::vector<std::string>& args);
  ConstantRightHandSideData(const int ncomp, const int comp, const double d);
  std::string GetName() const {return "constant";} 
  int GetNcomp() const {return _ncomp;}
  double operator()(int c, const Vertex2d& v)const;
  double operator()(int c, const Vertex3d& v)const;
};

/*-----------------------------------------*/

class OneComponentRightHandSideData : public DomainRightHandSide
{
 protected:
  
  int  ncomp, comp;  // ist die Komponente die Eins ist

 public:
  
  OneComponentRightHandSideData(int n, int c) : 
    DomainRightHandSide(), ncomp(n), comp(c) {}

  std::string GetName() const {return "one_onecomp";} 

  int GetNcomp() const { return ncomp;}

  double operator()(int c, const Vertex2d&)const 
    {
      if (c==comp) return 1.;
      return 0.;
    }
  double operator()(int c, const Vertex3d&)const 
    {
      if (c==comp) return 1.;
      return 0.;
    }
};

/*-----------------------------------------*/

class RectangleRightHandSideData : public DomainRightHandSide
{
  int      ncomp, comp;  // ist die Komponente die Eins ist
  double   x0, x1, y0, y1, z0, z1;

public:

  RectangleRightHandSideData(int n, int c, 
			     double xx0, double xx1, double yy0, double yy1) : 
    DomainRightHandSide(), ncomp(n), comp(c) 
    { 
      x0 = xx0; x1 = xx1; y0 = yy0; y1 = yy1;
      z0 = z1 = 0.;
    }
  RectangleRightHandSideData(int n, int c, 
			     double xx0, double xx1, double yy0, double yy1,
			     double zz0, double zz1) : 
    DomainRightHandSide(), ncomp(n), comp(c) 
    { 
      x0 = xx0; x1 = xx1; y0 = yy0; y1 = yy1;
      z0 = zz0; z1 = zz1;
    }

  std::string GetName() const {return "RectangleRightHandSideData";} 

  int GetNcomp() const {return ncomp;}

  double operator()(int c, const Vertex2d& V)const 
    {
      if (c!=comp) return 0.;
      if ((V.x()>x1) || (V.x()<x0)) return 0.;
      if ((V.y()>y1) || (V.y()<y0)) return 0.;
      return 1.;
    }
  double operator()(int c, const Vertex3d& V)const 
    {
      if (c!=comp) return 0.;
      if ((V.x()>x1) || (V.x()<x0)) return 0.;
      if ((V.y()>y1) || (V.y()<y0)) return 0.;
      if ((V.z()>z1) || (V.z()<z0)) return 0.;
      return 1.;
    }
};
}

#endif
