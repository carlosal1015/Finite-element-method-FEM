/**
*
* Copyright (C) 2004, 2007 by the Gascoigne 3D authors
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


#ifndef  __LOCAL1_h
#define  __LOCAL1_h

#include  "problemdescriptorbase.h"
#include  "laplace3d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidebyequation.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainmeanfunctional.h"

/*---------------------------------------------------*/

class PolynomialExactSolution : public Gascoigne::ExactSolution
{
  double quadratic(double x) const { return x*(1.-x);}

 public:
  
  std::string GetName()  const {return "PolynomialExactSolution";}
  int         GetNcomp() const { return 1; }
  double operator()(int c, const Gascoigne::Vertex3d& v) const
    {
      return quadratic(v.x()) * quadratic(v.y()) * quadratic(v.z());
    }
};

/*---------------------------------------------------*/

class ProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "Local";}
  void BasicInit(const Gascoigne::ParamFile* pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new Gascoigne::Laplace3d(GetParamFile());
    GetExactSolutionPointer() = new PolynomialExactSolution();
    GetRightHandSidePointer() = new Gascoigne::RightHandSideByEquation(GetEquation(), GetExactSolution());
    GetDirichletDataPointer() = new Gascoigne::DirichletDataByExactSolution(GetExactSolution());
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

class LocalDragFunctional : public virtual Gascoigne::ResidualFunctional
{
 public:
  LocalDragFunctional() : ResidualFunctional()
    {
      __comps.push_back(0);
      __cols.insert(9);
      __scales.push_back(1.);
      ExactValue() = 1./8.;

      __DD  = new Gascoigne::DirichletDataByColor(GetComps(),GetColors(),GetScales());
    }
  
  std::string GetName() const {
    return "LocalDrag";
  }
};

class LocalDomainFunctional : public virtual Gascoigne::AllDomainFunctional
{
 public:
  LocalDomainFunctional() : AllDomainFunctional(1,0)
    {
      ExactValue() = 0.02776989201546093;
    }

  std::string GetName() const {
    return "LocalDomain";
  }
};


#endif

