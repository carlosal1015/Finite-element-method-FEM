/**
*
* Copyright (C) 2008, 2009 by the Gascoigne 3D authors
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


#ifndef  __BenchSplitProblem_h
#define  __BenchSplitProblem_h

#include  "problemdescriptorbase.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include  "laplace2d.h"
#include  "domainrighthandside.h"
#include  "navierstokessplit.h"
#include  "filescanner.h"

/*-----------------------------------------*/

class BenchMarkSplitDirichletData : public Gascoigne::DirichletData
{
public:

  BenchMarkSplitDirichletData() {}
  std::string GetName() const {return "BenchMarkSplit";}
  void operator()(Gascoigne::DoubleVector& b, const Gascoigne::Vertex2d& v, int color) const 
  {
    b.zero();
/*     if (color!=80) */
/*       { */
/*         double high = 4.1; */
/* 	double vmax = 0.3; */
/* 	b[0] = vmax * Gascoigne::ParabelFunction(y,0.,high); */
/*       } */
  }
};

/*---------------------------------------------------*/

class VelocityRhs : public Gascoigne::DomainRightHandSide
{
  mutable const Gascoigne::TestFunction*  pressure;
  
  public  :

    VelocityRhs() : Gascoigne::DomainRightHandSide() {}
    std::string GetName() const { return "VelocityRhs" ;}
    int GetNcomp()        const { return 2;}
    void SetFemData(Gascoigne::FemData& d) const 
    { 
      std::map<const std::string,Gascoigne::FemFunction>::const_iterator it = d.find("pressure");
      pressure = &(it->second)[0];
    }
    double operator()(int c, const Gascoigne::Vertex2d& v)const 
    {
      if (c==0) return -pressure->x();
      else if (c==1) return -pressure->y();
      else abort();
    }
};

/*---------------------------------------------------*/

class VelocityProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "VelocityProblemDescriptor";}
  void BasicInit(const Gascoigne::ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new NavierStokesSplitLps2d(GetParamFile());
    GetRightHandSidePointer() = new VelocityRhs;
    GetDirichletDataPointer() = new BenchMarkSplitDirichletData;
    
    Gascoigne::ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

class PressureRhs : public Gascoigne::DomainRightHandSide
{
  mutable const Gascoigne::TestFunction  *v, *u;
  
  public  :

    PressureRhs() : Gascoigne::DomainRightHandSide() {}
    std::string GetName() const { return "PressureRhs" ;}
    int GetNcomp()        const { return 1;}
    void SetFemData(Gascoigne::FemData& d) const 
    { 
      std::map<const std::string,Gascoigne::FemFunction>::const_iterator it = d.find("velocity");
      v = &(it->second)[0];
      u = &(it->second)[1];
    }
    double operator()(int c, const Gascoigne::Vertex2d&)const 
    {
      return -( v->x() + u->y() ) / GetTimeStep();
    }
};

/*---------------------------------------------------*/

class NoDirichletBoundaries : public Gascoigne::BoundaryManager
{
 public:
  NoDirichletBoundaries() : Gascoigne::BoundaryManager() {}
 ~NoDirichletBoundaries() {}
 void BasicInit(const Gascoigne::ParamFile* pf) {}
};

/*---------------------------------------------------*/

class PressureProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "PressureProblemDescriptor";}
  void BasicInit(const Gascoigne::ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new Gascoigne::Laplace2d;
    GetRightHandSidePointer() = new PressureRhs;
    GetBoundaryManagerPointer() = new NoDirichletBoundaries;
    
    Gascoigne::ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

#endif
