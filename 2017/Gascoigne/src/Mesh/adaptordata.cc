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


#include "adaptordata.h"
#include "filescanner.h"

using namespace std;

/*****************************************************************/

namespace Gascoigne
{
void AdaptorData::ErrorMessage(const string& str1, const string& str2) const
{
  cout << "AdaptorData: " << str1 << " " << str2 << " ?" << endl;
  exit(1);
}

/*****************************************************************/

void AdaptorData::CheckEstimate() const
{
  if (_estimate=="none")      return;
  if (_estimate=="indicator") return;
  if (_estimate=="residuals") return;
  if (_estimate=="bench")     return;
  if (_estimate=="dual")      return;

  ErrorMessage("estimate",_estimate);
}

/*****************************************************************/

void AdaptorData::CheckAdaptor() const
{
  return;

//   if (_adaptor=="none")   return;
//   if (_adaptor=="nc")     return;
//   if (_adaptor=="opt")    return;
//   if (_adaptor=="diplomanden")   return;
//   if (_adaptor=="praktikum")     return;
//   if (_adaptor=="equal")  return;
//   if (_adaptor=="global") return;
//   if (_adaptor=="random") return;
//   if (_adaptor=="point")  return;
//   if (_adaptor=="bench3d") return;
//   if (_adaptor=="circle") return;
//   if (_adaptor=="pointpatch")  return;
//   if (_adaptor=="circlepatch") return;

//   ErrorMessage("adaptor",_adaptor);
}

/*****************************************************************/

void AdaptorData::CheckFunctional() const
{
  if (_functional=="none") return;
  if (_functional=="drag") return;
  if (_functional=="lift") return;
  if (_functional=="point") return;
  if (_functional=="mean") return;

  ErrorMessage("functional",_functional);
}

/*-------------------------------------------------*/

AdaptorData::AdaptorData() 
{
  mnodes    = 0;
  itol      = 0.01;
  icfactor  = 1.e8;
  irfactor  = 1.;
  idrfactor = 0.;
  idim      = 2;
  iorder    = 2;
  iglobal_conv = 1.;
  ilocal_conv = 2.;
}

/*-------------------------------------------------*/

void AdaptorData::read(const ParamFile* pf)
{
  DataFormatHandler DH;

  DH.insert("maxnodes"   ,& mnodes,10000000);
  DH.insert("rfactor"    ,& irfactor,1.);
  DH.insert("global_conv",& iglobal_conv,1.);
  DH.insert("local_conv" ,& ilocal_conv,2.);
  DH.insert("cfactor"    ,& icfactor,0.);
  DH.insert("estimate"   ,& _estimate,"none");
  DH.insert("adaptor"    ,& _adaptor ,"global");
  DH.insert("functional" ,& _functional,"none");
  DH.insert("randomp"    ,& _randomp,0.1);
  DH.insert("dimension"  ,& idim,2);

  FileScanner FS(DH,pf,"Adaptor");

//   CheckFunctional();
//   CheckAdaptor();
//   CheckEstimate();
}

/*-------------------------------------------------*/

void AdaptorData::reset()
{
  indr = inr = inc = 0;
  ieta = 0.; 
}
}
