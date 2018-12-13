/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#ifndef __baseq23dPatch_h
#define __baseq23dPatch_h

#include  "baseq23d.h"

#define NDOF   125
#define NDOF1d 5

namespace Gascoigne
{
/**************************************************/

  class BaseQ23dPatch : public Base3d
  {
    protected:
      bool         second;
      mutable int  pB;
      BaseQ23d     B2;
      nvector<int> perm;

      mutable nvector<double> dxx, dyy, dzz, dxy, dxz, dyz;

      void BasicInit()
      {
        N.resize(NDOF);
        DN.resize(NDOF);
        dxx.reservesize(NDOF);
        dyy.reservesize(NDOF);
        dzz.reservesize(NDOF);
        dxy.reservesize(NDOF);
        dxz.reservesize(NDOF);
        dyz.reservesize(NDOF);
      }

    public:
      BaseQ23dPatch() : Base3d(), second(false), perm(27)
      {
        BasicInit();

        perm[ 0] = 0;
        perm[ 1] = 1;
        perm[ 2] = 2;
        perm[ 3] = 5;
        perm[ 4] = 6;
        perm[ 5] = 7;
        perm[ 6] = 10;
        perm[ 7] = 11;
        perm[ 8] = 12;
        perm[ 9] = 25;
        perm[10] = 26;
        perm[11] = 27;
        perm[12] = 30;
        perm[13] = 31;
        perm[14] = 32;
        perm[15] = 35;
        perm[16] = 36;
        perm[17] = 37;
        perm[18] = 50;
        perm[19] = 51;
        perm[20] = 52;
        perm[21] = 55;
        perm[22] = 56;
        perm[23] = 57;
        perm[24] = 60;
        perm[25] = 61;
        perm[26] = 62;
      }

      int  n() const {return NDOF;}
      double phi   (int i) const {return N  [i];}
      double phi_x (int i) const {return DN [i].x();}
      double phi_y (int i) const {return DN [i].y();}
      double phi_z (int i) const {return DN [i].z();}
      double phi_xx(int i) const {assert(second); return dxx[i];}
      double phi_yy(int i) const {assert(second); return dyy[i];}
      double phi_zz(int i) const {assert(second); return dzz[i];}
      double phi_xy(int i) const {assert(second); return dxy[i];}
      double phi_xz(int i) const {assert(second); return dxz[i];}
      double phi_yz(int i) const {assert(second); return dyz[i];}
      const Vertex3d& phi_grad(int i) const {return DN [i];}

      void point(const Vertex3d& s) const
      {
        Vertex3d t(s);
        if( (s.x()<=0.5) && (s.y()<=0.5))
        {
          pB = 0;
        }
        else if( (s.x()>=0.5) && (s.y()<=0.5) )
        {
          pB = 2;
          t.x() -= 0.5;
        }
        else if( (s.x()<=0.5) && (s.y()>=0.5))
        {
          pB = 10;
          t.y() -= 0.5;
        }
        else
        {
          pB = 12;
          t.x() -= 0.5;
          t.y() -= 0.5;
        }
        if (s.z()>=0.5)
        {
          t.z() -= 0.5;
          pB += 50;
        }
        t *= 2.;

        B2.point(t);

        N.zero();
        for(int i=0; i<NDOF; i++)
        {
          DN[i].zero();
        }

        for (int i=0; i<27; i++)
        {
          int j = pB + perm[i];
          N [j] = B2.phi(i);
          DN[j] = B2.phi_grad(i);
          DN[j] *= 2.;
          if(second)
          {
            dxx[j] = 2.*B2.phi_xx(i);
            dyy[j] = 2.*B2.phi_yy(i);
            dzz[j] = 2.*B2.phi_zz(i);
            dxy[j] = 2.*B2.phi_xy(i);
            dxz[j] = 2.*B2.phi_xz(i);
            dyz[j] = 2.*B2.phi_yz(i);
          }
        }
      }
  };
}

#undef NDOF
#undef NDOF1d

#endif
