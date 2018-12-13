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


#ifndef  __TimeInfo_h
#define  __TimeInfo_h

#include  <string>

/*-----------------------------------------*/

namespace Gascoigne
{
class TimeInfo
{
protected:

  double _deltat, _time, _theta;
  double _tbegin, _tend;
  int    _iter, _neuler;
  double _ftscale[3], _fttheta[3];
  std::string _scheme, _actualscheme;
  double _Theta;

  int ftstep() const;

public:

  TimeInfo();

  void ReInit();
  void ReInit(double det);
  void BasicInit();

  double dt    () const;
  double theta () const;
  double oldrhs() const;
  double rhs   () const;

  void iteration(int i);
  void ReInit(double tb, double det, double te, const std::string& sch, int ne, double t);
  void ReInitTheta();
  void scale_timestep(double s) { _deltat *= s;}
  void stepback() { _time -= _deltat;}
  double time() const { return _time;}

  double time_begin() const { return _tbegin;}
  double time_end()   const { return _tend;}

  void ReInitBackward(int niter, double endtime);
  void iteration_backward(int i);
  void SpecifyScheme(int i);
  void RejectTimeStep(double d);
  void ScaleTimeStep(double d);
};
}

#endif
