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


#include  "timeinfo.h"

#include  <iostream>
#include  <stdio.h>
#include  <cassert>
#include  <math.h>

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
TimeInfo::TimeInfo()
{
  BasicInit();
}

/*-----------------------------------------*/

void TimeInfo::BasicInit()
{
  // fractional theta parameter
  //
  double gamma = 1.-sqrt(0.5);

  _ftscale[0] = gamma;
  _ftscale[1] = 1.-2.*gamma;
  _ftscale[2] = gamma;

  double alpha = (1.-2*gamma)/(1.-gamma);
  double beta = 1.-alpha;

  _fttheta[0] = alpha;
  _fttheta[1] = beta;
  _fttheta[2] = alpha;

  _theta = 1.;
  _scheme = "Euler";
  _actualscheme = "Euler";

  ReInit();
}

/*-----------------------------------------*/

void TimeInfo::ReInit()
{
  _time = 0.;
  _tbegin = _tend = _deltat = 0.;
  _iter = 0;
}

/*-----------------------------------------*/

void TimeInfo::ReInit(double det)
{
  _time   = _tbegin;
  _deltat = det;
  _iter = 0;
}

/*-----------------------------------------*/

int TimeInfo::ftstep() const
{
  int step = (_iter-_neuler)%3;
  while (step<0) step+=3;
  return step;
}

/*-----------------------------------------*/

void TimeInfo::ReInit(double tb, double te, double det, const string& sch, int ne, double tt)
{
  _tbegin = tb;
  _time   = _tbegin;
  _tend   = te;
  _deltat = det;
  _neuler = ne;
  _Theta  = tt;

  _scheme = sch;
  if (_neuler>0)
    _actualscheme = "Euler";
  else
    _actualscheme = sch;
	
  ReInitTheta();
}

/*-----------------------------------------*/

void TimeInfo::ReInitTheta()
{
  assert(_actualscheme=="CN" || _actualscheme=="Euler" || _actualscheme=="FractionalTheta" || _actualscheme=="Theta");

  if      (_actualscheme=="CN"   )           _theta = 0.5;
  else if (_actualscheme=="Euler")           _theta = 1.;
  else if (_actualscheme=="FractionalTheta") _theta = 0.;
  else if (_actualscheme=="Theta")           _theta = _Theta;
}

/*-----------------------------------------*/

double TimeInfo::theta() const 
{ 
  if (_actualscheme=="FractionalTheta") return _fttheta[ftstep()];
  return _theta;
}

/*-----------------------------------------*/

void TimeInfo::ScaleTimeStep(double d)
{
  _deltat *= d;
}

/*-----------------------------------------*/

double TimeInfo::dt() const 
{
  double h = 1.;
  if (_actualscheme=="FractionalTheta") h = _ftscale[ftstep()];
  
  assert((_deltat==0) || (_deltat>1.e-10));

  return _deltat * h;
}

/*-----------------------------------------*/

double TimeInfo::oldrhs() const 
{ 
  if (_actualscheme=="FractionalTheta")
    {
      if (ftstep()==1)
        {
          return 0.;
        }
      else
        {
          return 1./_fttheta[0];
        }
    }
  return (1.-_theta)/_theta;
}

/*-----------------------------------------*/

double TimeInfo::rhs() const 
{ 
  if (_actualscheme=="FractionalTheta")
    {
      if (ftstep()==2)
        {
          return 1./_fttheta[1];
        }
      else
        {
          return 0.;
        }
    }
  return 1.;
}

/*-----------------------------------------*/

void TimeInfo::ReInitBackward(int niter, double endtime)
{
  _iter = niter;
  _time = endtime;
}

/*-----------------------------------------*/

void TimeInfo::iteration_backward(int i)
{
  assert(i<=_iter);

  _iter = i;
  _time -= dt();

//  _actualscheme = _scheme;
//  if (_iter<=_neuler) _actualscheme = "Euler";

  cout << "\n============= " << _iter << " ========== " << _actualscheme;
  cout << " [t,dt] "<< _time << " " << dt() << "\t";
  cout << endl;
}

/*-----------------------------------------*/

void TimeInfo::RejectTimeStep(double d)
{
  _time -= dt();
  _deltat *= d;
}

/*-----------------------------------------*/

void TimeInfo::iteration(int i)
{
  assert(i>=_iter);

  _iter = i;
  _time += dt();

  cout << "\n============= " << _iter << " ========== " << _actualscheme;
  cout << " [t,dt] "<< _time << " " << dt() << "\t";
  cout << endl;
}

/*-----------------------------------------*/

void TimeInfo::SpecifyScheme(int i)
{
  if (i-1==_neuler)
    {
      cout << "Switching from Euler to " << _scheme << endl;
      _actualscheme = _scheme;
      ReInitTheta();
      if (_actualscheme=="FractionalTheta") _deltat *= 3.;
    }
}
}

/*-----------------------------------------*/

