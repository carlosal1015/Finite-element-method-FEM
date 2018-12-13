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


#include  "cfdblock3d.h"

using namespace std;

namespace Gascoigne
{
static CFDBlock3d  H;

/**********************************************************/

ostream& CFDBlock3d::print(ostream& os) const
{
  os << "CFD\n";
  os << s <<  " " << dx   << " " << dy << " " << dz << endl;
  os << gx << " " << laplx() << " " <<  0 <<  0 << endl;
  os << gy << " " <<  0   << " " << laply()<<  0 << endl;
  os << gz << " " <<  0   << 0 << " " << laplz()<<  endl;

  return os;
}

/**********************************************************/

void CFDBlock3d::zero()
{
  s=dx=dy=dz=gx=laplx()=laply()=laplz()=gy=gz=0.;
}

/**********************************************************/

void CFDBlock3d::DirichletRow (const vector<int>& cv)
{
  for (int i=0; i<cv.size(); i++) 
    {
      int c = cv[i];
      if      (c==1) gradx() = 0.;
      else if (c==2) grady() = 0.;
      else if (c==3) gradz() = 0.;
      if (c==1) laplx() = 0.;
      else if (c==2) laply() = 0.;
      else if (c==3) laplz() = 0.;
  }
}

/**********************************************************/

void CFDBlock3d::DirichletCol(const vector<int>& cv)
{
  for (int i=0; i<cv.size(); i++) 
    {
      int c = cv[i];
      if      (c==1) divx() = 0.;
      else if (c==2) divy() = 0.;
      else if (c==3) divz() = 0.;
      if (c==1) laplx() = 0.;
      else if (c==2) laply() = 0.;
      else if (c==3) laplz() = 0.;
  }
}

/**********************************************************/

void CFDBlock3d::DirichletDiag(const vector<int>& cv)
{
  for (int i=0; i<cv.size(); i++) 
    {
      int c = cv[i];
      if (c==1) laplx() = 1.;
      if (c==2) laply() = 1.;
      if (c==3) laplz() = 1.;
    }
}

/**********************************************************/

MatrixEntryType CFDBlock3d::operator()(int i,int j) const
{
  if (i==0)
    {
      assert(0<=j && j<=3);
      if (j==0) return s;
      if (j==1) return dx;
      if (j==2) return dy;
      if (j==3) return dz;
    }
  if (i==j) 
    {
      if (i==1) return laplx();
      if (i==2) return laply();
      if (i==3) return laplz();
    }
  if (j==0)
    {
      if (i==1) return gradx();
      if (i==2) return grady();
      if (i==3) return gradz();
    }
  return 0.;
}

/**********************************************************/

MatrixEntryType& CFDBlock3d::diag(int i)
{
  assert(0<=i && i<=3);
  if (i==0) return s;
  if (i==1) return laplx();
  if (i==2) return laply();
  if (i==3) return laplz();
  abort();
}

/**********************************************************/

void CFDBlock3d::adddiag(const nvector<double>& src, double l)
{
  s    += l * src[0];
  lapx += l * src[1];
  lapy += l * src[2];
  lapz += l * src[3];
}

/**********************************************************/

void CFDBlock3d::add(double d, const CFDBlock3d& A)
{
  s  += d * A.stab();
  dx += d * A.divx();
  dy += d * A.divy();
  dz += d * A.divz();
  gx += d * A.gradx();
  gy += d * A.grady();
  gz += d * A.gradz();
  lapx += d * A.laplx();
  lapy += d * A.laply();
  lapz += d * A.laplz();
}

/**********************************************************/

void CFDBlock3d::add(double d, const TimePattern& TP)
{
  s       += d*TP(0,0);
  laplx() += d*TP(1,1);
  laply() += d*TP(2,2);
  laplz() += d*TP(3,3);
}

/**********************************************************/

void CFDBlock3d::entry(const nmatrix<double>& E)
{
  s    += E(0,0);
  dx   += E(0,1);
  dy   += E(0,2);
  dz   += E(0,3);
  laplx() += E(1,1);
  laply() += E(2,2);
  laplz() += E(3,3);

  gx   += E(1,0);
  gy   += E(2,0);
  gz   += E(3,0);
}

/**********************************************************/

void CFDBlock3d::entry(int i, int j, const EntryMatrix& E, double d)
{
  s    += d*E(i,j,0,0);
  dx   += d*E(i,j,0,1);
  dy   += d*E(i,j,0,2);
  dz   += d*E(i,j,0,3);
  laplx() += d*E(i,j,1,1);
  laply() += d*E(i,j,2,2);
  laplz() += d*E(i,j,3,3);

  gx   += d*E(i,j,1,0);
  gy   += d*E(i,j,2,0);
  gz   += d*E(i,j,3,0);
}

/**********************************************************/

void CFDBlock3d::dual_entry(int i, int j, const EntryMatrix& E, double d)
{ 
  std::cerr << "\"CFDBlock3d::dual_entry\" not written!" << std::endl;
  abort();
//   s    += d*E(j,i,0,0);
//   dx   += d*E(j,i,1,0);
//   dy   += d*E(j,i,2,0);
//   dz   += d*E(j,i,3,0);
//   laplx() += d*E(i,j,1,1);
//   laply() += d*E(i,j,2,2);
//   laplz() += d*E(i,j,3,3);

//   gx   += d*E(j,i,0,1);
//   gy   += d*E(j,i,0,2);
//   gz   += d*E(j,i,0,3);
}

/**********************************************************/

void CFDBlock3d::operator *= (double d)
{
  s *= d;
  dx *= d; gx *= d;
  dy *= d; gy *= d;
  dz *= d; gz *= d;
  lapx *= d;
  lapy *= d;
  lapz *= d;
}

/**********************************************************/

void CFDBlock3d::operator *= (const CFDBlock3d& B)
{
  double news  = s*B.stab()  + dx*B.gradx() + dy*B.grady() + dz*B.gradz();
  double newdx = s*B.divx()  + dx*B.laplx(); 
  double newdy = s*B.divy()  + dy*B.laply(); 
  double newdz = s*B.divz()  + dz*B.laplz(); 

  double newgx = gradx()*B.stab() + laplx()*B.gradx();
  double newgy = grady()*B.stab() + laply()*B.grady();
  double newgz = gradz()*B.stab() + laplz()*B.gradz();

  double newlx = gradx()*B.divx();
  double newly = grady()*B.divy();
  double newlz = gradz()*B.divz();
  //double newl = (newlx+newly+newlz)/3. + laplx()*B.laplx(); 
  newlx += laplx()*B.laplx(); 
  newly += laply()*B.laply(); 
  newlz += laplz()*B.laplz(); 

  s  = news;
  dx = newdx;
  dy = newdy;
  dz = newdz;
  gx = newgx;
  gy = newgy;
  gz = newgz;
  laplx() = newlx;
  laply() = newly;
  laplz() = newlz;
}
/**********************************************************/

void CFDBlock3d::operator = (const CFDBlock3d& B)
{
  s  = B.stab();
  dx = B.divx();
  dy = B.divy();
  dz = B.divz();
  gx = B.gradx();
  gy = B.grady();
  gz = B.gradz();
  laplx() = B.laplx();
  laply() = B.laply();
  laplz() = B.laplz();
}

/**********************************************************/

void CFDBlock3d::operator += (const CFDBlock3d& B)
{
  s  += B.stab();
  dx += B.divx();
  dy += B.divy();
  dz += B.divz();
  gx += B.gradx();
  gy += B.grady();
  gz += B.gradz();
  laplx() += B.laplx();
  laply() += B.laply();
  laplz() += B.laplz();
}

/**********************************************************/

void CFDBlock3d::operator -= (const CFDBlock3d& B)
{
  s  -= B.stab();
  dx -= B.divx();
  dy -= B.divy();
  dz -= B.divz();
  gx -= B.gradx();
  gy -= B.grady();
  gz -= B.gradz();
  laplx() -= B.laplx();
  laply() -= B.laply();
  laplz() -= B.laplz();
}

/**********************************************************/

void CFDBlock3d::inverse()
{
  // A = ( C D )
  //     ( G L )
  //
  //  -1   ( C' D' )
  // A   = ( G' L' )
  //
  // mit S  = [ L-GC^{-1}D  ]^{-1} Schurkomplement
  //     L' =  S
  //     G' = -SGC^{-1}
  //     C' =  C^{-1}(I-DG')
  //     D' = -C^{-1}DS
  //

  double is = 1./s;
  double schurx = lapx - (gx*dx+gy*dy+gz*dz) * is;
  double schury = lapy - (gx*dx+gy*dy+gz*dz) * is;
  double schurz = lapz - (gx*dx+gy*dy+gz*dz) * is;
  schurx = 1./schurx;
  schury = 1./schury;
  schurz = 1./schurz;
  
  lapx = schurx;
  lapy = schury;
  lapz = schurz;

  double hx = -schurx * is;
  double hy = -schury * is;
  double hz = -schurz * is;
  gradx() *= hx;
  grady() *= hy;
  gradz() *= hz;

  s = (1.-dx*gradx()-dy*grady()-dz*gradz()) * is;

  dx *= hx;
  dy *= hy;
  dz *= hz;
}

/**********************************************************/

void CFDBlock3d::submult(const CFDBlock3d& B, const CFDBlock3d& C)
{
  // this -= B*C

  H  = B;
  H *= C;
  (*this) -= H;
}

/**********************************************************/

void CFDBlock3d::transpose()
{
  swap(dx,gx);
  swap(dy,gy);
  swap(dz,gz);
}

void CFDBlock3d::transpose(CFDBlock3d& A)
{ 
  swap(s , A.stab());
  swap(gx, A.divx());
  swap(gy, A.divy());
  swap(gz, A.divz());
  swap(dx, A.gradx());
  swap(dy, A.grady());
  swap(dz, A.gradz());
  swap(laplx(), A.laplx());
  swap(laply(), A.laply());
  swap(laplz(), A.laplz());
}


/**********************************************************/

void CFDBlock3d::vmult(iterator p) const
{
  double v0 = *p;
  double v1 = *(p+1);
  double v2 = *(p+2);
  double v3 = *(p+3);

  *p++ = s *v0+dx  *v1 + dy*v2 + dz*v3;
  *p++ = gradx()*v0+laplx()*v1;
  *p++ = grady()*v0+laply()*v2;
  *p   = gradz()*v0+laplz()*v3;
  
  p -=3;
}

/**********************************************************/

void CFDBlock3d::subtract(iterator p, const_iterator q0) const
{ 
  double a = s * *(q0) + dx * *(q0+1) + dy * *(q0+2) + dz * *(q0+3);
  *p++ -= a;
  
  a = gradx() * *q0  + laplx()* *(q0+1);
  *p++ -= a;
  
  a = grady() * *q0  + laply()* *(q0+2);
  *p++ -= a;

  a = gradz() * *q0  + laplz()* *(q0+3);
  *p -= a;

  p -= 3;
}
}
