/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#ifndef  __BoundaryRightHandSideByExactSolution_h
#define  __BoundaryRightHandSideByExactSolution_h


/////////////////////////////////////////////
////
////@brief
////  ... comments BoundaryRightHandSideByExactSolution

////
////
/////////////////////////////////////////////

#include  "boundaryrighthandside.h"
#include  "exactsolution.h"
#include  "equation.h"

class BoundaryRightHandSideByExactSolution : public Gascoigne::BoundaryRightHandSide
{
private:

  const Gascoigne::Equation*      _EQ;
  const Gascoigne::ExactSolution* _ES;

protected:


public:


//
////  Con(De)structor 
//
  
  BoundaryRightHandSideByExactSolution(const Gascoigne::Equation* eq, const Gascoigne::ExactSolution* es)
    : BoundaryRightHandSide(), _EQ(eq), _ES(es) { assert(es); assert(eq); }
  ~BoundaryRightHandSideByExactSolution() {}

  std::string GetName() const {return "BoundaryRightHandSideByExactSolution";}
  int GetNcomp() const { return _EQ->GetNcomp();}
  
  void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Gascoigne::Vertex2d& v, const Gascoigne::Vertex2d& n, int col) const{
    b[0] += ( _ES->x(0,v)*n.x()+_ES->y(0,v)*n.y() ) * N.m();
  }
  
 };


#endif
