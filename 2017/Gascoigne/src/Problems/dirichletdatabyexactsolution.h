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


#ifndef  __DirichletDataByExactSolution_h
#define  __DirichletDataByExactSolution_h

#include  "dirichletdata.h"
#include  "exactsolution.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class DirichletDataByExactSolution : public DirichletData
{
protected:

  const ExactSolution* ES;

public:

  DirichletDataByExactSolution(const ExactSolution* es) 
    : DirichletData(), ES(es)  {assert(es);}

  std::string GetName() const {return "ExactSolution";}
  
  void operator()(DoubleVector& b, const Vertex2d& v, int col)const{
    for(int c=0;c<b.size();c++) b[c] = (*ES)(c,v);
  }
  void operator()(DoubleVector& b, const Vertex3d& v, int col)const{
    for(int c=0;c<b.size();c++) b[c] = (*ES)(c,v);
  }


};
}

#endif
