/**
*
* Copyright (C) 2004, 2006, 2008 by the Gascoigne 3D authors
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


#ifndef __cginfo_h
#define __cginfo_h

#include  <string>

/*************************************************************/

namespace Gascoigne
{
class StatisticData
{
  double  _lastrate, _rate;
  int     _totaliter;

 public:
  StatisticData();

  double  rate     () const { return _rate; }
  double& rate     ()       { return _rate; }
  double& lastrate ()       { return _lastrate; }
  double  lastrate () const { return _lastrate; }
  int&    totaliter()       { return _totaliter; }
  int     totaliter() const { return _totaliter; }

  void reset();

  friend std::ostream& operator<<(std::ostream &s, const StatisticData& A);
};

/*************************************************************/

class ControlData
{
  std::string _status;
  int    _iteration;
  double _firstresidual, _residual, _aimedresidual;
  double _previousresidual;
  double _correction, _previouscorrection;

 public:

  ControlData();

  const std::string& status()     const { return _status;}
        std::string& status()           { return _status;}
  
  double& firstresidual()          { return _firstresidual;}
  double  firstresidual()    const { return _firstresidual;}
  double& residual     ()          { return _residual;}
  double  residual     ()    const { return _residual;}
  double& aimedresidual()          { return _aimedresidual;}
  double  aimedresidual()    const { return _aimedresidual;}
  double& previousresidual()       { return _previousresidual;}
  double  previousresidual() const { return _previousresidual;}
  int     iteration()        const { return _iteration;}
  int&    iteration()              { return _iteration;}
  double& correction     ()          { return _correction;}
  double  correction     ()    const { return _correction;}
  double& previouscorrection()       { return _previouscorrection;}
  double  previouscorrection() const { return _previouscorrection;}

  void reset();

  friend std::ostream& operator<<(std::ostream &s, const ControlData& A);
};

/*************************************************************/

class UserData
{
  std::string  _text;
  double _tol, _globaltol, _breaktol;
  int    _miniter,_maxiter, _printstep;

 public:

  double  tol()       const { return _tol;}
  double& tol()             { return _tol;}
  double  globaltol() const { return _globaltol;}
  double& globaltol()       { return _globaltol;}
  double  breaktol()  const { return _breaktol;}
  double& breaktol()        { return _breaktol;}
  int     miniter()   const { return _miniter;}
  int&    miniter()         { return _miniter;}
  int     maxiter()   const { return _maxiter;}
  int&    maxiter()         { return _maxiter;}
  int     printstep() const { return _printstep;}
  int&    printstep()       { return _printstep;}

  const std::string& text() const { return _text; }
        std::string& text()       { return _text; }

  friend std::ostream& operator<<(std::ostream &s, const UserData& A);
};

/*************************************************************/

class CGInfo
{
 protected:

  StatisticData SD;
  ControlData   CD;
  UserData      UD;
  
  void compute_reduction_rate();
  
 public:
  
  CGInfo(const std::string& txt);
  CGInfo(double f = 1.e-6, double r=1.e-14, int p = 10, 
	 int m = 100, const std::string& txt = "It: ");

  friend std::ostream& operator<<(std::ostream &s, const CGInfo& A);
  
  const StatisticData& statistics() const { return SD; }
        StatisticData& statistics()       { return SD; }
  const ControlData&   control()    const { return CD; }
        ControlData&   control()          { return CD; }
  const UserData&      user()       const { return UD; }
        UserData&      user()             { return UD; }
  
  bool check(double res, double corr=0.);
  void reset();
};
}

/*************************************************************/

#endif
