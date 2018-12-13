/**
*
* Copyright (C) 2006, 2011 by the Gascoigne 3D authors
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


#ifndef __baseq43d_h
#define __baseq43d_h

#include  <vector>
#include  <string>
#include  <utility>
#include  <cassert>
#include  "vertex.h"
#include  "numfixarray.h"
#include  "base3d.h"

#define NDOF   125
#define NDOF1d 5

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  Basis on the reference element

///
///
/////////////////////////////////////////////

class BaseQ43d : public Base3d
{

 protected:

  fixarray<NDOF1d,double>        a,b,c,d,e;

  double psi   (int i, double x) const { return a[i] + b[i]*x + c[i]*x*x + d[i]*x*x*x + e[i]*x*x*x*x;}
  double psi_x (int i, double x) const { return b[i] + 2.*c[i]*x + 3.*d[i]*x*x + 4.*e[i]*x*x*x;}
  double psi_xx(int i, double x) const { return 2.*c[i] + 6.*d[i]*x + 12.*e[i]*x*x;}

 public:

  BaseQ43d();

  int  n() const {return NDOF;}
  void point(const Vertex3d& s) const;

  double phi   (int i) const {return N  [i];}
  double phi_x (int i) const {return DN [i].x();}
  double phi_y (int i) const {return DN [i].y();}
  double phi_z (int i) const {return DN [i].z();}
  double phi_xx(int i) const {
    std::cerr << "\"BaseQ43d::phi_xx\" not written!" << std::endl;
    abort();
  }
  double phi_yy(int i) const {
    std::cerr << "\"BaseQ43d::phi_yy\" not written!" << std::endl;
    abort();
  }
  double phi_zz(int i) const {
    std::cerr << "\"BaseQ43d::phi_zz\" not written!" << std::endl;
    abort();
  }
  double phi_xy(int i) const {
    std::cerr << "\"BaseQ43d::phi_xy\" not written!" << std::endl;
    abort();
  }
  double phi_xz(int i) const {
    std::cerr << "\"BaseQ43d::phi_xz\" not written!" << std::endl;
    abort();
  }
  double phi_yz(int i) const {
    std::cerr << "\"BaseQ43d::phi_yz\" not written!" << std::endl;
    abort();
  }

  const Vertex3d&  phi_grad (int i) const {return DN [i];}
};
}

#undef NDOF
#undef NDOF1d

#endif
