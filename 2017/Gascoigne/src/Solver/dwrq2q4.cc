/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#include "dwrq2q4.h"
#include "dwrfemq2.h"
#include "piq2.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

DwrQ2Q4::DwrQ2Q4(SolverInterface &S) : _S(S)
{
  _P = _S.GetProblemDescriptor();
  _D = _S.GetDiscretization();
}

/**********************************************************/

DiscretizationInterface* DwrQ2Q4::CreateOtherDiscretization() const
{
  DiscretizationInterface* D;

  if(_S.GetMesh()->dimension()==2)
  {
    D = new DwrFemQ2Q42d;
  }
  else
  {
    D = new DwrFemQ2Q43d;
  }

  D->BasicInit(_S.GetParamfile());
  return D;
}

/**********************************************************/

double DwrQ2Q4::ScalarProduct(DoubleVector &eta, const GlobalVector &f, const GlobalVector &z) const
{
  for(int i=0; i<z.n(); i++)
  {
    for(int c=0; c<z.ncomp(); c++)
    {
      eta[i] += fabs(f(i,c)*z(i,c));
    }
  }
  return z * f;
}

/**********************************************************/

double DwrQ2Q4::ScalarProduct(DoubleVector &eta, const VectorInterface &gf, const VectorInterface &gz) const
{
  const GlobalVector& f = _S.GetGV(gf);
  const GlobalVector& z = _S.GetGV(gz);

  return ScalarProduct(eta,f,z);
}

/**********************************************************/

double DwrQ2Q4::ScalarProductWithFluctuations(DoubleVector& eta, const VectorInterface &gf, const VectorInterface &gz) const
{
  const GlobalVector& f = _S.GetGV(gf);

  GlobalVector dz(f.ncomp(),f.n());
  PiQ2 pi;
  pi.Init(_S.GetMesh());
  pi.vmult(dz,_S.GetGV(gz));

  return ScalarProduct(eta,f,dz);
}

/**********************************************************/

void DwrQ2Q4::PrimalResidualsHigher(VectorInterface &gf, const VectorInterface &gu)
{
  _S.GetGV(gf).zero();

  // only necessary if z has additional Dirichlet bc compared to u
  //
  _S.Rhs (gf,  -0.5);
  _S.Form(gf,gu,0.5);

  DiscretizationInterface* D = CreateOtherDiscretization();

  _S.SetDiscretization(*D,true);

  _S.Rhs (gf,    0.5);
  _S.Form(gf,gu,-0.5);

  _S.SetDiscretization(*_D);
  delete D;
}

/**********************************************************/

void DwrQ2Q4::DualResidualsHigher(VectorInterface &gf, const VectorInterface &gu, const VectorInterface &gz, const ProblemDescriptorInterface &PDI)
{
  _S.GetGV(gf).zero();
  // dual problem
  _S.SetProblem(PDI);
  _S.AddNodeVector("u",gu);

  // standard residual
  //
  {
    _S.Rhs (gf,  -0.5);
    _S.Form(gf,gz,0.5);
  }
  // residual respect Q4 test functions
  //
  {
    DiscretizationInterface* D = CreateOtherDiscretization();
    _S.SetDiscretization(*D,true);

    _S.Rhs (gf,    0.5);
    _S.Form(gf,gz,-0.5);

    _S.SetDiscretization(*_D);
    delete D;
  }

  _S.DeleteNodeVector("u");
  _S.SetProblem(*_P);
}

/**********************************************************/

double DwrQ2Q4::Estimator(DoubleVector &eta, VectorInterface &gf, const VectorInterface &gu, const VectorInterface &gz, const ProblemDescriptorInterface &PDI)
{
  double rho=0, rhostern=0;
  eta.resize(_S.GetGV(gz).n());

  DualResidualsHigher(gf,gu,gz,PDI);
  rhostern = ScalarProductWithFluctuations(eta,gf,gu);

  PrimalResidualsHigher(gf,gu);
  rho      = ScalarProductWithFluctuations(eta,gf,gz);

  return rho + rhostern;
}

/**********************************************************/
}

