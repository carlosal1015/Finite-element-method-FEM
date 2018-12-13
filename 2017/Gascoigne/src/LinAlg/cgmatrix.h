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


#ifndef __cgmatrix_h
#define __cgmatrix_h


namespace Gascoigne
{
template <class MT,class VT,class MEM, class INFO>
inline int CgMatrix(const MT& M, VT& x, const VT& b, MEM& mem, INFO& info)
{
  
  VT& g  = mem[0];
  VT& h  = mem[1];
  VT& d  = mem[2];
  VT& Ad = mem[3];

  int    it=0,reached = 0;
  double gh,alpha,beta;
  double     res;

  M.vmulteq(g,x);
  g.sequ(1.,-1.,b);
  res = g.norm();

  if(res==0.)      return 0;

  d = g;
  //M.precondition(d,g);

  gh  =  g*d;
  gh *=  -1.;

  for(it=0;!reached;it++)
  {
    M.vmulteq(Ad,d);

    alpha = d*Ad;
    alpha = gh/alpha;

    g.add(alpha,Ad);
    x.add(alpha,d );
    res = g.norm();

    reached = info.check_residual(it,res);
    if (reached) continue;
    
    h = g;
    //M.precondition(h,g);
    
    beta = gh;
    gh   = g*h;
    beta = gh/beta;
    
    d.sequ(beta,-1.,h);
  }
  if (reached<0) return 1;
  return 0;
}
}

#endif
