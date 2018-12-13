/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#ifndef  __DomainMeanFunctional_h
#define  __DomainMeanFunctional_h

#include  "domainfunctional.h"
#include  <set>

/*-----------------------------------------*/

namespace Gascoigne
{
class AllDomainFunctional  : public virtual DomainFunctional
{
protected:

  int  _comp, _ncomp;

public:

  AllDomainFunctional(int nc, int c) { _ncomp = nc; _comp = c; }
  ~AllDomainFunctional() {}

  std::string GetName() const {return "AllDomainFunctional";}

  int    GetNcomp() const {return _ncomp;}
  int    GetComp()  const {return _comp;}

  double J(const FemFunction& U, const Vertex2d& v) const;
  double J(const FemFunction& U, const Vertex3d& v) const;
};

/*-----------------------------------------*/

class SubDomainFunctional  : public AllDomainFunctional
{
protected:

  double  _x0, _x1, _y0, _y1, _z0, _z1;

public:

  SubDomainFunctional(int nc, int c) : AllDomainFunctional(nc,c) {};
  ~SubDomainFunctional() {}

  std::string GetName() const {return "SubDomainFunctional";}

  void SetCoordinates(double x0, double x1, double y0, double y1, double z0=0, double z1=0);

  double J(const FemFunction& U, const Vertex2d& v) const;
  double J(const FemFunction& U, const Vertex3d& v) const;
};
}

#endif
