/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#include "stdmultilevelsolverdata.h"
#include "filescanner.h"
#include "stringutil.h"

using namespace std;

/**********************************************************/

namespace Gascoigne
{
StdMultiLevelSolverData::~StdMultiLevelSolverData()
{
}

/**********************************************************/

void StdMultiLevelSolverData::BasicInit(const ParamFile *param)
{
  _countresidual = 0; 
  
  double prec_tol, prec_globaltol;
  int    prec_maxiter, prec_pstep;

  DataFormatHandler DFH;
  DFH.insert("linearsolve",         &_linearsolve,        "mg");
  DFH.insert("nonlinearsolve",      &_nonlinearsolve,     "newton");
  DFH.insert("solver",              &_solver,             "stat");
  DFH.insert("mgomega",             &_mgomega,            1.);
  DFH.insert("coarselevel",         &_coarselevel,        0);
  DFH.insert("mgtype",              &_mgtype,             "V");

  DFH.insert("prec_tol",      &prec_tol,      1.e-12);
  DFH.insert("prec_globaltol",&prec_globaltol,1.e-12);
  DFH.insert("prec_maxiter",  &prec_maxiter,  1);
  DFH.insert("prec_pstep",    &prec_pstep,    0);
  DFH.insert("gmressize",     &_gmresmemsize, 10);

  DFH.insert("save_nonlinear_comp_residuals",  &_i_save_nonlinear_comp_residuals, 0);
  DFH.insert("save_linear_comp_residuals",     &_i_save_linear_comp_residuals,    0);
  DFH.insert("show_nonlinear_comp_residuals",  &_i_show_nonlinear_comp_residuals, 0);
  DFH.insert("show_linear_comp_residuals",     &_i_show_linear_comp_residuals,    0);
  DFH.insert("show_comp_residual_names",       &_i_show_comp_residual_names,      0);

  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(param,"MultiLevelSolver");

  if ((_mgtype!="V") && (_mgtype!="W") && (_mgtype!="F"))
  {
    _mgtype = "V";
  }
  precinfo.user().tol()       = prec_tol;
  precinfo.user().globaltol() = prec_globaltol;
  precinfo.user().maxiter()   = prec_maxiter;
  precinfo.user().printstep() = prec_pstep;
  precinfo.user().text()      = "PrecInfo";
}
}
