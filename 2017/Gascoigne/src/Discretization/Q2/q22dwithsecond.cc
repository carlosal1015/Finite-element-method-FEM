/**
*
* Copyright (C) 2004, 2005, 2008 by the Gascoigne 3D authors
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


#include "q22dwithsecond.h"
#include "integratorwithsecond.h"

using namespace std;

namespace Gascoigne
{

/**********************************************************/

void Q22dWithSecond::BasicInit(const ParamFile* paramfile)
{
  if (GetIntegrator()==NULL)
    GetIntegratorPointer() =  new IntegratorWithSecond<2>;
  assert(GetIntegrator());

  typedef FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond> FEWithSecond;

  if (GetFem()==NULL)
    GetFemPointer() =  new FEWithSecond;
  assert(GetFem());

  Q22d::BasicInit(paramfile);
}

/* ----------------------------------------- */

double Gascoigne::Q22dWithSecond::EstimateSecond(DoubleVector& eta, const GlobalVector& u, double dd) const
{
  DoubleVector d(9);
  d[0] = 1.; d[1] = 2.; d[2] = 1.; d[3] = 2.; d[4] = 4.; d[5] = 2.; d[6] = 1.; d[7] = 2.; d[8] = 1.;

  DoubleVector eeta(GetMesh()->nnodes(),0.);

  nmatrix<double> T;

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocal(__U,u,iq);
    dynamic_cast<const IntegratorWithSecond<2>*>(GetIntegrator())->EstimateSecond(__F,*GetFem(),__U);

    IntVector indices = GetLocalIndices(iq);

    for(int ii=0; ii<indices.size(); ii++) 
    {
      int i = indices[ii];
      for(int c=0;c<__F.ncomp(); c++)
      {
        eeta[i] += d[ii]*sqrt(__F(0,c))*dd;
      }
    }
  }

  double sum=0;
  for(int i=0; i<GetMesh()->nnodes(); i++)
  {
    sum += eeta[i] * eeta[i];
  }

  eta.add(1,eeta);

  return sum;
}
}
