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


#ifndef __baseq12dPatch_h
#define __baseq12dPatch_h

#include  "baseq12d.h"

#define NDOF   9
#define NDOF1d 3

namespace Gascoigne
{
/**************************************************/

class BaseQ12dPatch : public Base2d
{
 protected:

  mutable int        pB;
  BaseQ12d           B1;

  mutable Vertex2d   bn, bt;
  nvector<int>       perm;
  mutable nvector<double> dxy;

  void BasicInit() {
      N.resize(NDOF);  
      DN.resize(NDOF);  
/*       dxy.resize(NDOF);   */
    }

 public:
  
  BaseQ12dPatch() : perm(4) 
  {
    BasicInit();
    
    perm[0] = 0;
    perm[1] = 1;
    perm[2] = 3;
    perm[3] = 4;
  }

  int  n() const {return NDOF;}
  double phi   (int i) const {return N  [i];}
  double phi_x (int i) const {return DN [i].x();}
  double phi_y (int i) const {return DN [i].y();}
  double phi_xx(int i) const {return 0.;}
  double phi_xy(int i) const {
    std::cerr << "\"BaseQ12dPatch::phi_xy\" not written!" << std::endl;
    abort();
//     return dxy[i];
  }
  double phi_yy(int i) const {return 0.;}
  const Vertex2d&  phi_grad (int i) const {return DN [i];}

  void point(const Vertex2d& s) const
    {
      Vertex2d t(s);
      if( (s.x()<=0.5) && (s.y()<=0.5) )
	{
	  pB = 0;
	}
      else if( (s.x()>=0.5) && (s.y()<=0.5) )
	{
	  pB = 1;
	  t.x() -= 0.5;
	}
      else if( (s.x()<=0.5) && (s.y()>=0.5) )
	{
	  pB = 2;
	  t.y() -= 0.5;
	}
      else
	{
	  pB = 3;
	  t.x() -= 0.5;
	  t.y() -= 0.5;
	}
      t *= 2.;

      B1.point(t);

      N.zero();
      for(int i=0;i<NDOF;i++)  DN[i].zero();

      int offset = perm[pB];
      for(int i=0; i<4; i++)
	{
	  int j = offset + perm[i];
	  N[j]  = B1.phi(i);
	  DN[j] = B1.phi_grad(i);  
	  DN[j] *= 2.;
/* 	  dxy[j] = B1.phi_xy(i);   */
/* 	  dxy[j] *= 2.; */
	}
    }
};
}

#undef NDOF
#undef NDOF1d

#endif
