/**
*
* Copyright (C) 2008 by the Gascoigne 3D authors
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


#ifndef  __BenchMarkProblem_h
#define  __BenchMarkProblem_h

#include  "problemdescriptorbase.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include  "navierstokeslps2d.h"

/*-----------------------------------------*/

class BenchMarkDirichletData : public Gascoigne::DirichletData
{
public:

  BenchMarkDirichletData() {}
  std::string GetName() const {return "Bench";}
  void operator()(Gascoigne::DoubleVector& b, const Gascoigne::Vertex2d& v, int color) const 
  {
    double y = v.y();

    b.zero();
    if (color!=80)
      {
        double high = 4.1;
	double vmax = 0.3;
	b[1] = vmax * Gascoigne::ParabelFunction(y,0.,high);
      }
  }
};

/*---------------------------------------------------*/

class BenchMarkProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "NavierStokesBenchmark";}
  void BasicInit(const Gascoigne::ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new Gascoigne::NavierStokesLps2d(GetParamFile());
    GetDirichletDataPointer() = new BenchMarkDirichletData;
    
    Gascoigne::ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

#endif
