/**
*
* Copyright (C) 2004, 2005, 2006, 2008 by the Gascoigne 3D authors
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


#include "solverinfos.h"
#include <cassert>
#include "filescanner.h"

using namespace std;

namespace Gascoigne
{

SolverInfos::~SolverInfos()
{
  map<string,CGInfo *>::iterator p = __L.begin();
  while(p != __L.end())
  { 
    if(p->second) 
    {
      delete p->second; 
      p->second = NULL;
    } 
    p++;
  }
  
  map<string,NLInfo *>::iterator q = _NL.begin();
  while(q != _NL.end())
  { 
    if(q->second) 
    {
      delete q->second; 
      q->second = NULL;
    } 
    q++;
  }
}

/*---------------------------------------------------------------*/

CGInfo& SolverInfos::GetLInfo(string s) const
{
  map<string,CGInfo *>::const_iterator iter = __L.find(s);
  assert(iter!=__L.end());
  return *iter->second; 
}

/*---------------------------------------------------------------*/

NLInfo& SolverInfos::GetNLInfo(string s) const
{
  map<string,NLInfo *>::const_iterator iter = _NL.find(s);
  assert(iter!=_NL.end());
  return *iter->second; 
}

/*---------------------------------------------------------------*/

void SolverInfos::BasicInit(const ParamFile *param)
{
  __L["State"]   = new CGInfo();
  _NL["State"]   = new NLInfo(GetLInfo());
  
  DataFormatHandler DFH;

  double linear_tol, linear_globaltol;
  int    linear_miniter, linear_maxiter, linear_pstep;

  double nonlinear_tol, nonlinear_globaltol, nonlinear_rho, nonlinear_increase;
  int    nonlinear_miniter, nonlinear_maxiter, nonlinear_pstep, nonlinear_damp;

  DFH.insert("linearsolve",         &_linearsolve,        "mg");
  DFH.insert("linear_tol",          &linear_tol,          1.e-2);
  DFH.insert("linear_globaltol",    &linear_globaltol,    1.e-12);
  DFH.insert("linear_miniter",      &linear_miniter,       0);
  DFH.insert("linear_maxiter",      &linear_maxiter,      10);
  DFH.insert("linear_pstep",        &linear_pstep,        0);

  DFH.insert("nonlinear_tol",       &nonlinear_tol,       1.e-4);
  DFH.insert("nonlinear_globaltol", &nonlinear_globaltol, 1.e-12);
  DFH.insert("nonlinear_rho",       &nonlinear_rho,       0.001);
  DFH.insert("nonlinear_miniter",   &nonlinear_miniter,    0);
  DFH.insert("nonlinear_maxiter",   &nonlinear_maxiter,   10);
  DFH.insert("nonlinear_pstep",     &nonlinear_pstep,     0);
  DFH.insert("nonlinear_damp",      &nonlinear_damp,      4);
  DFH.insert("nonlinear_increase",  &nonlinear_increase,  1.e3);

  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(param,"MultiLevelSolver");

  GetLInfo().user().tol()                = linear_tol;
  GetLInfo().user().globaltol()          = linear_globaltol;
  GetLInfo().user().miniter  ()          = linear_miniter;
  GetLInfo().user().maxiter  ()          = linear_maxiter;
  GetLInfo().user().printstep()          = linear_pstep;

  GetLInfo().user().text() = "Lin:";
                                         
  GetNLInfo().user().tol()               = nonlinear_tol;
  GetNLInfo().user().globaltol()         = nonlinear_globaltol;
  GetNLInfo().user().rho()               = nonlinear_rho;
  GetNLInfo().user().miniter()           = nonlinear_miniter;
  GetNLInfo().user().maxiter()           = nonlinear_maxiter;
  GetNLInfo().user().printstep()         = nonlinear_pstep;
  GetNLInfo().user().maxrelax()          = nonlinear_damp;
  GetNLInfo().user().maxresincrease()    = nonlinear_increase;
}

/*---------------------------------------------------------------*/
}

