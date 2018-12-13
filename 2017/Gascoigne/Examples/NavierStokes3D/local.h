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


#ifndef  __local_h
#define  __local_h

#include  "stdloop.h"
#include  "navierstokeslps3d.h"
#include  "stokeslps3d.h"
#include  "dirichletdata.h"
#include  "meshagent.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include  "problemdescriptorbase.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

class BenchMarkDirichletData : public DirichletData
{
protected:
  double vmax;
public:
  BenchMarkDirichletData() {
    vmax = 0.45;
  }
  std::string GetName() const {return "Bench";}
  void operator()(DoubleVector& b, const Vertex3d& v, int color) const {

    double x = v.x();  double y = v.y(); double z = v.z();
    
    b.zero();
    if (x==0.)
      {
	double high = 4.1;
	b[1] = vmax * ParabelFunction(y,0.,high) * ParabelFunction(z,0.,high);
      }
  }
};

/* ----------------------------------------- */

class ProblemDescriptor : public ProblemDescriptorBase
{
public:
    
    std::string GetName() const {return "Local";}
    void BasicInit(const ParamFile* pf) {
      GetParamFilePointer() = pf;
      GetEquationPointer() = new StokesLps3d(GetParamFile());
      GetDirichletDataPointer() = new BenchMarkDirichletData();
      ProblemDescriptorBase::BasicInit(pf);
    }
};


/* ----------------------------------------- */

class RunderKreis : public BoundaryFunction<2>
{
  double squareradius;
  Vertex2d center;
public :

  std::string GetName() const { return "RunderKreis";}

  void BasicInit(Vertex2d c, double r) 
    {
      center = c; 
      squareradius = r;
    }
  double operator()(const Vertex2d& c) const 
    {
      double r = - squareradius;
      for (int i=0; i<2; i++)
	{
	  double dx = c[i]-center[i];
	  r += dx * dx;
	}
      return r;
    }
};

/*---------------------------------------------------*/

class BenchMarkMeshAgent : public MeshAgent
{
 protected:
  
  RunderKreis RK;

 public:
  
  BenchMarkMeshAgent() : MeshAgent()
    {
      double r = 0.25;
      Vertex2d v(2.,2.);
      RK.BasicInit(v,r);
      AddShape(80,&RK);
    }
};

/* ----------------------------------------- */

class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC) 
    {
      GetMeshAgentPointer() = new BenchMarkMeshAgent;
      StdLoop::BasicInit(paramfile, PC);
    }
};


#endif
