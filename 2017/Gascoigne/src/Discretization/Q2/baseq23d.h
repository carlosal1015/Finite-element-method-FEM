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


#ifndef __baseq23d_h
#define __baseq23d_h

#include  <vector>
#include  <string>
#include  <utility>
#include  <cassert>
#include  "vertex.h"
#include  "numfixarray.h"
#include  "base3d.h"

#define NDOF   27
#define NDOF1d 3

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  Basis on the reference element

///
///
/////////////////////////////////////////////

class BaseQ23d : public Base3d
{

 protected:

  fixarray<NDOF1d,double>        a,b,c;
  
  double psi   (int i, double x) const { return a[i] + b[i]*x + c[i]*x*x;}
  double psi_x (int i, double x) const { return b[i] + 2.*c[i]*x;       }
  double psi_xx(int i, double x) const { return 2.*c[i];       }

 public:
  
  BaseQ23d();

  int  n() const {return NDOF;}
  void point(const Vertex3d& s) const;

  double phi   (int i) const {return N  [i];}
  double phi_x (int i) const {return DN [i].x();}
  double phi_y (int i) const {return DN [i].y();}
  double phi_z (int i) const {return DN [i].z();}
  double phi_xx(int i) const {
    std::cerr << "\"BaseQ23d::phi_xx\" not written!" << std::endl;
    abort();
  }
  double phi_yy(int i) const {
    std::cerr << "\"BaseQ23d::phi_yy\" not written!" << std::endl;
    abort();
  }
  double phi_zz(int i) const {
    std::cerr << "\"BaseQ23d::phi_zz\" not written!" << std::endl;
    abort();
  }
  double phi_xy(int i) const {
    std::cerr << "\"BaseQ23d::phi_xy\" not written!" << std::endl;
    abort();
  }
  double phi_xz(int i) const {
    std::cerr << "\"BaseQ23d::phi_xz\" not written!" << std::endl;
    abort();
  }
  double phi_yz(int i) const {
    std::cerr << "\"BaseQ23d::phi_yz\" not written!" << std::endl;
    abort();
  }

  const Vertex3d&  phi_grad (int i) const {return DN [i];}
};
}

#undef NDOF
#undef NDOF1d

#endif
