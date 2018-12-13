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


#ifndef __baseq12d_h
#define __baseq12d_h

#include  <vector>
#include  <string>
#include  <utility>
#include  "vertex.h"
#include  "numfixarray.h"
#include  "base2d.h"

#define NDOF   4
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

class BaseQ12d : public Base2d
{
 protected:

  fixarray<2,double>          a,b;
  mutable DoubleVector dxy;

  void BasicInit();

  double psi_x(int i, double x) const { return b[i]; }

 public:
  
  BaseQ12d();

  double psi  (int i, double x) const { return a[i] + b[i]*x;}
  int    n()                    const { return NDOF;}
  void   point(const Vertex2d& s) const;

  double phi   (int i) const {return N  [i];}
  double phi_x (int i) const {return DN [i].x();}
  double phi_y (int i) const {return DN [i].y();}
  double phi_xx(int i) const {return 0.;}
  double phi_yy(int i) const {return 0.;}
/*   double phi_xy(int i) const {return 0.;} */
  double phi_xy(int i) const
  {
    std::cerr << "\"BaseQ12d::phi_xy\" not written!" << std::endl;
    abort();
//     return dxy[i];
  }

  const Vertex2d &  phi_grad (int i) const {return DN [i];}
};
}

#undef NDOF
#undef NDOF1d

#endif
