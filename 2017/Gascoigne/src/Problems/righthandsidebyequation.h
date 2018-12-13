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


#ifndef  __RightHandSideByEquation_h
#define  __RightHandSideByEquation_h

#include  "domainrighthandside.h"
#include  "exactsolution.h"
#include  "equation.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class RightHandSideByEquation : public DomainRightHandSide
{
protected:

  const Equation*      _EQ;
  const ExactSolution* _ES;

public:

  RightHandSideByEquation(const Equation* eq, const ExactSolution* es)
    : DomainRightHandSide(), _EQ(eq), _ES(es) { assert(es); assert(eq);}
  
  std::string GetName() const { return "RightHandSideByEquation";} 
  int GetNcomp() const { return _EQ->GetNcomp();}

  double operator()(int c, const Vertex2d& v)const 
    {
      int n = _EQ->GetNcomp();
      DoubleVector b(n,0.);
      FemFunction U(n);
      for (int i=0; i<n; i++)
	{
	  U[i].m() = (*_ES)(i,v);
	  U[i].x() = _ES->x(i,v);
	  U[i].y() = _ES->y(i,v);
	  U[i].D() = _ES->xx(i,v)+_ES->yy(i,v);
	}
      _EQ->OperatorStrong(b,U);
      if (GetTimeStep()>0.)
	{
	  double eps = 1.e-6;
	  DoubleVector ut(n,0.);
	  double time = GetTime();
	 _ES->SetTime(time+0.5*eps,GetTimeStep());
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] = (*_ES)(i,v);
	    }
	  _ES->SetTime(time-0.5*eps,GetTimeStep());
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] -= (*_ES)(i,v);
	    }
	  _ES->SetTime(time,GetTimeStep());
	  b.add(1./eps,ut);
	}
      return b[c];
    }
  double operator()(int c, const Vertex3d& v)const 
    {
      int n = _EQ->GetNcomp();
      DoubleVector b(n,0.);
      FemFunction U(n);
      for (int i=0; i<n; i++)
	{
	  U[i].m() = (*_ES)(i,v);
	  U[i].x() = _ES->x(i,v);
	  U[i].y() = _ES->y(i,v);
	  U[i].z() = _ES->z(i,v);
	  U[i].D() = _ES->xx(i,v)+_ES->yy(i,v)+_ES->zz(i,v);
	}
      _EQ->OperatorStrong(b,U);
      if (GetTimeStep()>0.)
	{
	  double eps = 1.e-6;
	  DoubleVector ut(n,0.);
	  double time = GetTime();
	 _ES->SetTime(time+0.5*eps,GetTimeStep());
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] = (*_ES)(i,v);
	    }
	  _ES->SetTime(time-0.5*eps,GetTimeStep());
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] -= (*_ES)(i,v);
	    }
	  _ES->SetTime(time,GetTimeStep());
	  b.add(1./eps,ut);
	}
      return b[c];
    }

};
}

#endif
