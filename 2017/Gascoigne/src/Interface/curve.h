/**
*
* Copyright (C) 2007 by the Gascoigne 3D authors
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


#ifndef  __Curve_h
#define  __Curve_h

#include "vertex.h"
#include "nvector.h"


namespace Gascoigne
{
/*-----------------------------------------*/


class Curve
{
 private:

 protected:

  
 public:

  Curve() {}
  virtual ~Curve() {}

  virtual int GetNcomp() const=0;

  virtual double X(double t) const = 0;
  virtual double Y(double t) const = 0;

  virtual double DX(double t) const = 0;
  virtual double DY(double t) const = 0;
  
  virtual Vertex2d operator()(double t) const
    {
      Vertex2d x;
      x.x()=X(t);
      x.y()=Y(t);
      return x;
    }
  
  virtual double NormD(double t) const {return sqrt(DX(t)*DX(t) + DY(t)*DY(t));}

  virtual void Normal(nvector<double>& n,double t) const
    {
      n.resize(2);
      double d = NormD(t);
      assert(d>=0);
      n[0] =  DY(t)/d;
      n[0] = -DX(t)/d;
    }

  virtual void Vertices(nvector<Vertex2d>& V,const nvector<double>& t) const
    {
      assert(V.size()==t.size());
      for(int i=0;i<t.size();++i)
	{
	  V[i].x() = X(t[i]);
	  V[i].y() = Y(t[i]);
	}
    }
};

}

#endif
