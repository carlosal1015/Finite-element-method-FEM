/**
*
* Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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

#include  "domainrighthandside.h"
#include  "equation.h"
#include  "filescanner.h"
#include  "problemdescriptorbase.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

class LocalEquation : public Equation
{
 private:
  mutable double _visc;
  mutable double _us, _vs;
  mutable double _h, _k, _r;
 public:
  LocalEquation(const ParamFile* paramfile) {
    DataFormatHandler DFH;
    DFH.insert("visc",&_visc,1.);
    DFH.insert("h",&_h,0.4);
    DFH.insert("k",&_k,2.);
    DFH.insert("r",&_r,0.3);
    FileScanner FS(DFH,paramfile,"Equation");
    
    _us = (_r*_h)/(1.-_r);
    _vs = (1.-_us)*(_h+_us);
    cerr << "****** us, vs = **** " << _us << " " << _vs << endl; 
  }

  double GetUs() const {return _us;}
  double GetVs() const {return _vs;}

  std::string GetName() const { return "Local";}
  int  GetNcomp      () const {return 2;}
  void SetTimePattern(TimePattern& P) const {
    P.reservesize(GetNcomp(),GetNcomp(),0.);
    P(0,0) = 1.;
    P(1,1) = 1.;
  }

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const {
    b[0] += _visc* (U[0].x()*N.x()+U[0].y()*N.y());
    b[1] += _visc* (U[1].x()*N.x()+U[1].y()*N.y());
    
    double u = U[0].m();
    double v = U[1].m();
    double s = u/(u+_h);
    
    b[0] += N.m() * (-u*(1.-u) + s * v);
    b[1] += N.m() * (_k*_r* v - _k* s * v);
  }

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const {
    A(0,0) += _visc* (M.x()*N.x()+M.y()*N.y());
    A(1,1) += _visc* (M.x()*N.x()+M.y()*N.y());
    
    double u = U[0].m();
    double v = U[1].m();
    double MM = M.m()*N.m();
    
    double s = u/(u+_h);
    double t = _h/((u+_h)*(u+_h));
    
    A(0,0) += ( (-1.+2.*u) + t*v ) * MM;
    A(0,1) += s * MM;
    A(1,0) += -_k * t*v * MM;
    A(1,1) += ( _k*_r - _k*s ) * MM;
  }
};

/* ----------------------------------------- */

class LocalInitialCondition : public DomainRightHandSide
{
private:
  mutable double _us, _vs;
public:
  LocalInitialCondition(const LocalEquation* EQ) {
    _us = EQ->GetUs();
    _vs = EQ->GetVs();
  }
  
  std::string GetName() const { return "Local";}
  int GetNcomp() const {return 2;}  
  double operator()(int c, const Vertex2d& v) const {
  double x = v.x();
  double y = v.y();
  double eps1 = 2.e-7;
  double eps2 = 3.e-5;
  double eps3 = 1.2e-4;
  if(c==0)
    {
      double dist = - eps1*(x-0.1*y-225)*(x-0.1*y-675);
      return _us + dist;
    }
  else if(c==1)
    {
      return _vs - eps2*(x-450)-eps3*(y-150);
    }
  abort();
  }
};

/* ----------------------------------------- */

class ProblemDescriptor : public ProblemDescriptorBase
{
public:
    
    std::string GetName() const {return "Local";}
    void BasicInit(const ParamFile* pf) {
      GetParamFilePointer() = pf;
      GetEquationPointer() = new LocalEquation(GetParamFile());
      const LocalEquation* LEQ = dynamic_cast<const LocalEquation*>(GetEquation());
      GetInitialConditionPointer() = new LocalInitialCondition(LEQ);

      ProblemDescriptorBase::BasicInit(pf);
    }
};

#endif
