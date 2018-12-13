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


#ifndef  __HeatProblem_h
#define  __HeatProblem_h

#include  "problemdescriptorbase.h"
#include  "boundaryfunction.h"

/*-----------------------------------------*/

class HeatDirichletData : public Gascoigne::DirichletData
{
public:

  HeatDirichletData() {}

  std::string GetName() const {return "Heat";}

  void operator()(Gascoigne::DoubleVector& b, const Gascoigne::Vertex2d& v, int color) const {

    b.zero();
    if (color!=80)
      {
        b[0] = 10;
        b[1] = 20;
      }
  }
};

/* ----------------------------------------- */

class HeatEquation : public Gascoigne::Equation
{
 private:

  mutable double _visc;
  mutable double _us, _vs;
  mutable double _h, _k, _r;

 public:

  HeatEquation(const Gascoigne::ParamFile* paramfile) 
    {
      Gascoigne::DataFormatHandler DFH;
      DFH.insert("visc",&_visc,1.);
      DFH.insert("h",&_h,0.4);
      DFH.insert("k",&_k,2.);
      DFH.insert("r",&_r,0.3);
      Gascoigne::FileScanner FS(DFH,paramfile,"Equation");
      
      _us = (_r*_h)/(1.-_r);
      _vs = (1.-_us)*(_h+_us);
      //std::cout << "****** us, vs = **** " << _us << " " << _vs << std::endl;
    }

  std::string GetName() const { return "HeatEquation";}
  int  GetNcomp      () const {return 2;}
  void SetTimePattern(Gascoigne::TimePattern& P) const {
    P.reservesize(GetNcomp(),GetNcomp(),0.);
    P(0,0) = 1.;
    P(1,1) = 1.;
  }

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
	    const Gascoigne::TestFunction& N) const 
  {
    b[0] += _visc* (U[0].x()*N.x()+U[0].y()*N.y());
    b[1] += _visc* (U[1].x()*N.x()+U[1].y()*N.y());

    double u = U[0].m();
    double v = U[1].m();
    double s = u/(u+_h);

    b[0] += N.m() * (-u*(1.-u) + s * v);
    b[1] += N.m() * (_k*_r* v - _k* s * v);
  }

  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
	      const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const 
  {
    A(0,0) += _visc* (M.x()*N.x()+M.y()*N.y());
    A(1,1) += _visc* (M.x()*N.x()+M.y()*N.y());

    double u = U[0].m();
    double v = U[1].m();
    double MM = M.m()*N.m();

    double s = u/(u+_h);
    double t = _h/((u+_h)*(u+_h));

    A(0,0) += MM * ( (-1.+2.*u) + t*v );
    A(0,1) += MM * s;
    A(1,0) -= MM * _k * t*v;
    A(1,1) += MM * _k *(_r-s);
  }
};

/*---------------------------------------------------*/

class HeatProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "HeatProblemDescriptor";}
  void BasicInit(const Gascoigne::ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new HeatEquation(GetParamFile());
    GetDirichletDataPointer() = new HeatDirichletData;
    
    Gascoigne::ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

#endif
