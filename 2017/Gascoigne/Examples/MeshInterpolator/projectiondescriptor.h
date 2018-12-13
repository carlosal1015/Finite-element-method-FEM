/**
*
* Copyright (C) 2005 by the Gascoigne 3D authors
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


#ifndef  __ProjectionDescriptor_h
#define  __ProjectionDescriptor_h

#include  "problemdescriptorbase.h"

/*---------------------------------------------------*/

class ProjectionEquation : public Gascoigne::Equation
{
 public:

  ProjectionEquation() : Gascoigne::Equation() { }
  ~ProjectionEquation() { }

  int GetNcomp() const { return 1; }
  std::string GetName() const { return "ProjectionEquation"; }

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
	    const Gascoigne::TestFunction& N) const
  {
    b[0] += U[0].m() * N.m();
  }
  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
	      const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const
  {
    A(0,0) += M.m() * N.m();
  }
};

/*---------------------------------------------------*/

class ProjectionProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "Projection";}
  void BasicInit(const Gascoigne::ParamFile* pf) 
  {
    GetEquationPointer() = new ProjectionEquation();
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

#endif
