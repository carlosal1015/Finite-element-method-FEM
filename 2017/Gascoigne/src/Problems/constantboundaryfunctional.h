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


#ifndef  __ConstantBoundaryFunctional_h
#define  __ConstantBoundaryFunctional_h


#include  "boundaryfunctional.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class ConstantBoundaryFunctional : public BoundaryFunctional
{
protected:

  int              comp;
  std::set<int>    colors;
  double           value;

public:


  ConstantBoundaryFunctional();
  ConstantBoundaryFunctional(const std::vector<std::string>& args);
  ~ConstantBoundaryFunctional();
  void Construct(const std::vector<std::string>& args);
  
  std::string GetName() const {return "ConstantBoundaryFunctional";}

  std::set<int> GetColors() const {return colors;}

  void AddColor(int    c) {colors.insert(c);}
  void SetComp (int    c) {comp =c;}
  void SetValue(double v) {value=v;}

  double J(const FemFunction& U, const Vertex2d& v) const;

};
}

#endif
