/**
*
* Copyright (C) 2004, 2005, 2011 by the Gascoigne 3D authors
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


#ifndef __baseq13d_h
#define __baseq13d_h

#include  "base3d.h"

#define NDOF   8
#define NDOF1d 2

/**************************************************/

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  Basis on the reference element

///
///
/////////////////////////////////////////////

class BaseQ13d : public Base3d
{
 protected:

  typedef  std::pair<int,int> pint;

  fixarray<2,double>          a,b;

  mutable DoubleVector dxy, dxz, dyz, dxyz;

  void BasicInit()
    {
      N.resize(NDOF);  DN.resize(NDOF);  
      a[0] = 1.;  b[0] = -1.;
      a[1] = 0.;  b[1] =  1.;
/*       dxy.reservesize(NDOF); */
/*       dxz.reservesize(NDOF); */
/*       dyz.reservesize(NDOF); */
/*       dxyz.reservesize(NDOF); */
    }
  double psi_x(int i, double x)      const { return b[i]; }

 public:
  
  BaseQ13d();

  double psi  (int i, double x)      const { return a[i] + b[i]*x;}
  int    n()                         const { return NDOF;}
  double phi   (int i)               const { return N [i];}
  double phi_x (int i)               const { return DN[i].x();}
  double phi_y (int i)               const { return DN[i].y();}
  double phi_z (int i)               const { return DN[i].z();}
  double phi_xx(int i) const {return 0.;}
  double phi_yy(int i) const {return 0.;}
  double phi_zz(int i) const {return 0.;}
  double phi_xy(int i) const {
    std::cerr << "\"BaseQ13d::phi_xy\" not written!" << std::endl;
    abort();
//     return dxy[i];
  }
  double phi_xz(int i) const {
    std::cerr << "\"BaseQ13d::phi_xz\" not written!" << std::endl;
    abort();
//     return dxz[i];
  }
  double phi_yz(int i) const {
    std::cerr << "\"BaseQ13d::phi_yz\" not written!" << std::endl;
    abort();
//     return dyz[i];
  }
  double phi_xyz(int i) const {
    std::cerr << "\"BaseQ13d::phi_xyz\" not written!" << std::endl;
    abort();
//     return dxyz[i];
  }
  const Vertex3d &  phi_grad (int i) const { return DN[i];}

  void point(const Vertex3d& s) const;
};
}

#undef NDOF
#undef NDOF1d

#endif
