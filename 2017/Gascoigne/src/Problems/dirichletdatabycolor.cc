/**
*
* Copyright (C) 2004, 2007, 2008 by the Gascoigne 3D authors
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


#include  "dirichletdatabycolor.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  DirichletDataByColor::DirichletDataByColor(nvector<int> comps, const set<int>& cl, nvector<double> s)
    : DirichletData(), __cols(cl), __comps(comps), __scales(s)
  {
    assert(__comps.size()==__scales.size());
    assert(__comps.size()>0);
    assert(__cols.size()>0);
  }

  DirichletDataByColor::DirichletDataByColor(int comps, std::set<int>& cl, double s) : __cols(cl)
  {
    __comps .push_back(comps);
    __scales.push_back(s);
  }

  /*-----------------------------------------*/

  DirichletDataByColor::DirichletDataByColor(const vector<string>& args)
  {
    bool ok=true;
    
    int n = args.size(); 
    if (n<5) ok=false;
    int ncol = atoi(args[0].c_str());
    if (n<4+ncol) ok = false;
    for (int i=0;i<ncol;++i)
      __cols.insert(atoi(args[i+1].c_str()));
    int ncomp = atoi(args[ncol+1].c_str());
    if (n!=2*ncomp+ncol+2) ok = false;
    for (int i=0;i<ncol;++i)
      {
	__comps.push_back( atoi(args[2*i+ncol+2].c_str()));
	__scales.push_back(atof(args[2*i+ncol+3].c_str()));
      }
    
    if (!ok)
      {
	cerr << "DirichletDataByColor::DirichletDataByColor: Usage" << endl
	     << "\t ncols_[cols]_ncomps_[comp_weight] " << endl;
	abort();
      }
  }

  /*-----------------------------------------*/

  void DirichletDataByColor::operator()
    (DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();

    if(__cols.find(color)!=__cols.end()) 
      for (int i=0;i<__comps.size();++i)
	b[__comps[i]] = __scales[i];
  }

  /*-----------------------------------------*/

  void DirichletDataByColor::operator()
    (DoubleVector& b, const Vertex3d& v, int color) const
  {
    b.zero();
    if(__cols.find(color)!=__cols.end()) 
      for (int i=0;i<__comps.size();++i)
	b[__comps[i]] = __scales[i];
  }

  /*-----------------------------------------*/

  set<int> DirichletDataByColor::preferred_colors()const 
  {
    return __cols;
  }
}
