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


#ifndef __base2d_h
#define __base2d_h

#include  "base.h"
#include  "gascoigne.h"

/**************************************************/

namespace Gascoigne
{
class Base2d : public Base
{
 protected:

  mutable DoubleVector                N;
  mutable std::vector<Vertex2d>       DN;

  mutable Vertex2d  bn, bt;

 public:
  
  const Vertex2d*  normal2d() const { return &bn;}
  const Vertex2d*  tangent2d() const { return &bt;}

  void point_boundary(int ie, const Vertex1d& s1) const
    {
      Vertex2d s;
      if     (ie==0)      
	{
	  s.x() = s1.x(); s.y() = 0.;     
	  bn.x() =  0.; bn.y() = -1.;
	  bt.x() =  1.; bt.y() =  0.;
	}
      else if(ie==1)      
	{
	  s.x() = 1.    ; s.y() = s1.x(); 
	  bn.x() =  1.; bn.y() =  0.;
	  bt.x() =  0.; bt.y() =  1.;
	}
      else if(ie==2)      
	{
	  s.x() = s1.x(); s.y() = 1.;     
	  bn.x() =  0.; bn.y() =  1.;
	  bt.x() = -1.; bt.y() =  0.;
	}
      else                
	{
	  s.x() = 0.    ; s.y() = s1.x(); 
	  bn.x() = -1.; bn.y() =  0.;
	  bt.x() =  0.; bt.y() = -1.;
	}
      point(s);
    }
  
};
}

#endif
