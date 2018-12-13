/**
*
* Copyright (C) 2004, 2008 by the Gascoigne 3D authors
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


#include  "ncadaptor.h"
#include  "compareclass.h"
#include  "giota.h"
#include  "filescanner.h"
#include  <limits>


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
NCAdaptor::NCAdaptor(const ParamFile* paramfile, const DoubleVector& _eta) 
  :  eta(_eta)
{
  DataFormatHandler DH;

#ifdef __OLDCOMPILER__
  DH.insert("n"   ,& _n,INT_MAX);
#else
  DH.insert("n"   ,& _n,numeric_limits<int>::max());
#endif

  DH.insert("p"  ,& _p,0.9);
  FileScanner FS(DH, paramfile, "Adaptor");

  etasum = accumulate(eta.begin(),eta.end(),0.);
}

/*-----------------------------------------*/

void NCAdaptor::refine(IntVector& ref, IntVector& coars) const
{
  int n = eta.size();
  IntVector C(n); 
  iota(C.begin(),C.end(),0);
  typedef CompareObjectBigToSmall<DoubleVector >  CoC;
  sort(C.begin(),C.end(),CoC(eta));

  double eta_limit = _p*etasum;

  int    i    = 0;
  double done = 0.;

  ref.clear();
  while((done<eta_limit) && (i<_n))
    {
      done += fabs(eta[C[i]]);
      ref.push_back(C[i]);
      i++;
    }
}
}
