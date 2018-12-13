/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#include "meshagent.h"
#include "stdsolver.h"
#include "stdmultilevelsolver.h"
#include "gascoigne.h"
#include "paramfile.h"
#include "q12d.h"

using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

class LocalSolver : public StdSolver
{
  public:

    LocalSolver() : StdSolver() {};
    ~LocalSolver() {}

    double Integral(const VectorInterface& gu) const
    {
      const GlobalVector& u = GetGV(gu);
      HNAverage(gu);
      assert(GetDiscretization()->HNZeroCheck(u)==0);

      nvector<double> dst = _PF.IntegrateVector(u);
      HNZero(gu);
      return dst[0];
    }
    DiscretizationInterface* NewDiscretization(int dimension, const string& discname)
    {
      return new Q12d;
    }
    void NewMesh(int l, const MeshInterface* MP)
    {
      StdSolver::NewMesh(l,MP);
      GetDiscretization()->InitFilter(_PF);
    }
    void BasicInit(int level, const ParamFile* paramfile, const int dimension)
    {
      StdSolver::BasicInit(level,paramfile,dimension);
    }
};

/*---------------------------------------------------*/

double f(double x, double t)
{
  double xx = 0.5 + 0.25*cos(2*M_PI*t);
  double yy = 0.5 + 0.25*sin(2*M_PI*t);
  double u = 1./(1.+(x-xx)*(x-xx)+yy*yy);
  double val = u*u*yy;
  return val;
}

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  int niter = 10;

  MeshAgent M;
  M.BasicInit(&paramfile);

  VectorInterface gu("u");

  for (int iter=1; iter<=niter; iter++)
  {
    LocalSolver S;

    const MeshInterface* MI = M.GetMesh(0);
    S.BasicInit(iter,&paramfile,MI->dimension());
    S.NewMesh(iter,MI);

    S.ReInitVector(gu,1);
    GlobalVector& u=S.GetGV(gu);
    u.resize(MI->nnodes());

    for (int i=0; i<u.n(); i++)
    {
      double x = M.GetMesh()->vertex2d(i).x();
      double y = M.GetMesh()->vertex2d(i).y();
      u(i,0) = f(x,y);
    }

    double val = S.Integral(gu);
    cout.precision(16);
    cout << iter << "\t" << u.n() << "\t" << val << endl;

    M.global_refine(1);
  }
}
