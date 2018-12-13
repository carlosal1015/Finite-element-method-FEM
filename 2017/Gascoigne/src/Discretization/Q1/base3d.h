/**
*
* Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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


#ifndef __base3d_h
#define __base3d_h

#include  "base.h"
#include  "gascoigne.h"

/**************************************************/

namespace Gascoigne
{
class Base3d : public Base
{
 protected:

  mutable DoubleVector             N;
  mutable std::vector<Vertex3d>       DN;
  mutable fixarray<2,int>             face;
  mutable Vertex3d  bn, bt;

 public:
  
  Base3d()  {}
  const Vertex3d*  normal3d() const {return &bn;}
  const Vertex3d*  tangent3d() const {
    std::cerr << "\"Base3d::tangent3d\" not written!" << std::endl;
    abort();
  }
  const fixarray<2,int>* faces() const {return &face;}
  
  void point_boundary(int ie, const Vertex2d& s1) const
    {
      Vertex3d s;
      if     (ie==0)      
	{
	   s.x() = s1.x(); s.y() = s1.y(); s.z() = 0.;
	  bn.x() = 0.;    bn.y() = 0.;    bn.z() = -1.;
	  face[0] = 0; face[1] = 1;
	}
      else if(ie==1)      
	{
	   s.x() = 1.;  s.y() = s1.y(); s.z() = s1.x();
	  bn.x() = 1.; bn.y() = 0.;    bn.z() = 0.;
	  face[0] = 2; face[1] = 1;
	}
      else if(ie==2)      
	{
	   s.x() = s1.x(); s.y() = 1.;  s.z() = s1.y();
	  bn.x() = 0.;    bn.y() = 1.; bn.z() = 0.;
	  face[0] = 0; face[1] = 2;
	}
      else if(ie==3)
	{
	   s.x() = 0.; s.y() = s1.y(); s.z() = 1.-s1.x();
	  bn.x() = -1.; bn.y() =  0.; bn.z() = 0.;
	  face[0] = 1; face[1] = 2;
	}
      else if(ie==4)      
	{
	   s.x() = s1.x(); s.y() = 0.;   s.z() = 1.-s1.y();
	  bn.x() = 0.;    bn.y() = -1.; bn.z() = 0.;
	  face[0] = 0; face[1] = 2;
	}
      else
	{
	   s.x() = 1.-s1.x(); s.y() = s1.y(); s.z() = 1.;
	  bn.x() = 0.;       bn.y() = 0.;    bn.z() = 1.;
	  face[0] = 0; face[1] = 1;
	}
      point(s);
    }
};
}

#endif
