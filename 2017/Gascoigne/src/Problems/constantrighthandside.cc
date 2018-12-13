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


#include  "constantrighthandside.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

ConstantRightHandSideData::ConstantRightHandSideData
(const int ncomp, const int comp, const double d) 
  : DomainRightHandSide()
{
  _ncomp = ncomp;
  _comp  = comp;
  _d     = d;
}

ConstantRightHandSideData::ConstantRightHandSideData
(const vector<string>& args) 
  : DomainRightHandSide()
{
  if(args.size()!=3)
    {
      cerr << "ConstantRightHandSideData::ConstantRightHandSideData()\n";
      cerr << "wrong number of args (3): " << args.size() << endl;
      abort();
    }
  _ncomp = atoi(args[0].c_str());
  _comp  = atoi(args[1].c_str());
  _d     = atof(args[2].c_str());
}

/*-----------------------------------------*/

double ConstantRightHandSideData::operator()(int c, const Vertex2d& v)const 
{
  if(c==_comp) return _d;
  return 0.;
}

/*-----------------------------------------*/

double ConstantRightHandSideData::operator()(int c, const Vertex3d& v)const 
{
  if(c==_comp) return _d;
  return 0.;
}
}
