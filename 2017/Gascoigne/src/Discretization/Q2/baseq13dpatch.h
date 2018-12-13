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


#ifndef __baseq13dPatch_h
#define __baseq13dPatch_h

#include  "baseq13d.h"

#define NDOF  27
#define NDOF1d 3

namespace Gascoigne
{
/**************************************************/

class BaseQ13dPatch : public Base3d
{
 protected:

  mutable int                pB;
  BaseQ13d           B1;
  nvector<int>       perm;

/*   mutable nvector<double>             N; */
/*   mutable std::vector<Vertex3d>   DN; */

  mutable nvector<double> dxy, dxz, dyz, dxyz;

  void BasicInit() {
      N.resize(NDOF);  
      DN.resize(NDOF);  
/*       dxy.reservesize(NDOF); */
/*       dxz.reservesize(NDOF); */
/*       dyz.reservesize(NDOF); */
/*       dxyz.reservesize(NDOF); */
    }

 public:
  
  BaseQ13dPatch() : perm(8)        
  {
    BasicInit();
    perm[0] = 0;
    perm[1] = 1;
    perm[2] = 3;
    perm[3] = 4;
    perm[4] = 9;
    perm[5] = 10;
    perm[6] = 12;
    perm[7] = 13;
  }

  int  n() const {return NDOF;}
  double phi   (int i) const {return N  [i];}
  double phi_x (int i) const {return DN [i].x();}
  double phi_y (int i) const {return DN [i].y();}
  double phi_z (int i) const {return DN [i].z();}
  double phi_xx(int i) const {return 0.;}
  double phi_yy(int i) const {return 0.;}
  double phi_zz(int i) const {return 0.;}
  double phi_xy(int i) const {
    std::cerr << "\"BaseQ13dPatch::phi_xy\" not written!" << std::endl;
    abort();
//     return dxy[i];
  }
  double phi_xz(int i) const {
    std::cerr << "\"BaseQ13dPatch::phi_xz\" not written!" << std::endl;
    abort();
//     return dxz[i];
  }
  double phi_yz(int i) const {
    std::cerr << "\"BaseQ13dPatch::phi_yz\" not written!" << std::endl;
    abort();
//     return dyz[i];
  }
  double phi_xyz(int i) const {
    std::cerr << "\"BaseQ13dPatch::phi_xyz\" not written!" << std::endl;
    abort();
//     return dxyz[i];
  }
  const Vertex3d&  phi_grad (int i) const {return DN [i];}

  void point(const Vertex3d& s) const
    {
      Vertex3d t(s);
      if( (s.x()<=0.5) && (s.y()<=0.5))
	{
	  pB = 0;
	}
      else if( (s.x()>=0.5) && (s.y()<=0.5) )
	{
	  pB = 1;
	  t.x() -= 0.5;
	}
      else if( (s.x()<=0.5) && (s.y()>=0.5))
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
      if (s.z()>=0.5)
	{
	  t.z() -= 0.5;
	  pB += 4;
	}
      t *= 2.;

      B1.point(t);

      N.zero();
      for(int i=0;i<NDOF;i++)  DN[i].zero();

      int offset = perm[pB];

      for (int i=0; i<8; i++)
	{
	  int j = perm[i] +offset;
	  N [j] = B1.phi(i);
	  DN[j] = B1.phi_grad(i);
	  DN[j] *= 2.;
/* 	  dxy[j] = 2.*B1.phi_xy(i); */
/* 	  dxz[j] = 2.*B1.phi_xz(i); */
/* 	  dyz[j] = 2.*B1.phi_yz(i); */
/* 	  dxyz[j] = 2.*B1.phi_xyz(i); */
	}

/*       int offset = 0; */
/*       if      (pB==1) offset = 1; */
/*       else if (pB==2) offset = 3; */
/*       else if (pB==3) offset = 4; */
/*       else if (pB==4) offset = 9; */
/*       else if (pB==5) offset = 10; */
/*       else if (pB==6) offset = 12; */
/*       else if (pB==7) offset = 13; */
      
/*       N[ 0+offset] = B1.phi(0); */
/*       N[ 1+offset] = B1.phi(1); */
/*       N[ 4+offset] = B1.phi(2); */
/*       N[ 3+offset] = B1.phi(3); */
/*       N[ 9+offset] = B1.phi(4); */
/*       N[10+offset] = B1.phi(5); */
/*       N[13+offset] = B1.phi(6); */
/*       N[12+offset] = B1.phi(7); */
      
/*       DN[ 0+offset] = B1.phi_grad(0); */
/*       DN[ 1+offset] = B1.phi_grad(1); */
/*       DN[ 4+offset] = B1.phi_grad(2); */
/*       DN[ 3+offset] = B1.phi_grad(3); */
/*       DN[ 9+offset] = B1.phi_grad(4); */
/*       DN[10+offset] = B1.phi_grad(5); */
/*       DN[13+offset] = B1.phi_grad(6); */
/*       DN[12+offset] = B1.phi_grad(7); */

/*       for(int i=0;i<NDOF;i++)  DN[i] *= 2.; */
    }
};
}

#undef NDOF
#undef NDOF1d

#endif
