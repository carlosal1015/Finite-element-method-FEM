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


#ifndef __baseq22dPatch_h
#define __baseq22dPatch_h

#include  "baseq22d.h"

#define NDOF   25
#define NDOF1d 5

namespace Gascoigne
{
/**************************************************/

  class BaseQ22dPatch : public Base2d
  {
    protected:
      bool                 second;
      mutable int          pB;
      BaseQ22d             B2;
      nvector<int>         perm;
      mutable DoubleVector dxx,dxy,dyy;

      void BasicInit()
      {
        N.resize(NDOF);
        DN.resize(NDOF);
        dxy.resize(NDOF);
      }

    public:
      BaseQ22dPatch() : Base2d(), second(false), perm(9)
      {
        BasicInit();

        perm[0] = 0;
        perm[1] = 1;
        perm[2] = 2;
        perm[3] = 5;
        perm[4] = 6;
        perm[5] = 7;
        perm[6] = 10;
        perm[7] = 11;
        perm[8] = 12;
      }

      int  n() const {return NDOF;}
      double phi   (int i) const {return N  [i];}
      double phi_x (int i) const {return DN [i].x();}
      double phi_y (int i) const {return DN [i].y();}
      double phi_xx(int i) const {assert(second); return dxx[i];}
      double phi_xy(int i) const {assert(second); return dxy[i];}
      double phi_yy(int i) const {assert(second); return dyy[i];}
      const Vertex2d& phi_grad(int i) const {return DN [i];}

      void point(const Vertex2d& s) const
      {
        Vertex2d t(s);
        if( (s.x()<=0.5) && (s.y()<=0.5) )
        {
          pB = 0;
        }
        else if( (s.x()>=0.5) && (s.y()<=0.5) )
        {
          pB = 2;
          t.x() -= 0.5;
        }
        else if( (s.x()<=0.5) && (s.y()>=0.5) )
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
        t *= 2.;

        B2.point(t);

        N.zero();
        for(int i=0; i<NDOF; i++)
        {
          DN[i].zero();
        }

        for(int i=0; i<9; i++)
        {
          int j  = pB + perm[i];
          N[j]   = B2.phi(i);
          DN[j]  = B2.phi_grad(i);
          DN[j] *= 2.;
          if(second)
          {
            dxx[j] = 2.*B2.phi_xx(i);
            dxy[j] = 2.*B2.phi_xy(i);
            dyy[j] = 2.*B2.phi_yy(i);
          }
        }
      }
  };
}

#undef NDOF
#undef NDOF1d

#endif
