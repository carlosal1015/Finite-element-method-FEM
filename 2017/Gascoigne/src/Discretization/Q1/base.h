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


#ifndef __base_h
#define __base_h

#include  "vertex.h"
#include  <string>

/**************************************************/

namespace Gascoigne
{
class Base
{
 protected:
  
  void error(const std::string& s) const
    {
      std::cout << "Base::" << s << " not written !\n"; 
      abort();
    }

 public:

  virtual ~Base(){}
  
  virtual const Vertex2d*  normal2d() const {return NULL;}
  virtual const Vertex2d*  tangent2d() const {return NULL;}

  virtual const Vertex3d*  normal3d() const {return NULL;}
  virtual const Vertex3d*  tangent3d() const {return NULL;}
  virtual const fixarray<2,int>* faces() const
    { error("faces");return NULL;}

  virtual double psi(int i, double x) const
    { error("psi"); return 0;}

  virtual void point(const Vertex2d&) const
    { error("point");}

  virtual void point(const Vertex3d&) const
    { error("point");}

  virtual void point_boundary_2d(int, const Vertex1d&) const
    { error("point_boundary");}

  virtual void point_boundary_3d(int, const Vertex2d&) const
    { error("point_boundary");}  
};
}

#endif
