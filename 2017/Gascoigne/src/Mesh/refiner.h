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


#ifndef __refiner_h
#define __refiner_h

#include "meshinterface.h"
#include "nmatrix.h"

/*************************************************************/

namespace Gascoigne
{
class PointRefiner
{
  const MeshInterface& H;
  const Vertex2d& V;

  nmatrix<double> A;
  bool VertexInQuad(int);

 public:

  PointRefiner(const MeshInterface& h, const Vertex2d& v) 
    : H(h), V(v), A(2,2) {}

  void BuildCellList(std::vector<int>&);
};

/*************************************************************/

class CircleRefiner 
{
  const MeshInterface&  H;
  const Vertex3d&  V;
  double R;

  bool QuadOnRadius(int) const;

 public:

  CircleRefiner(const MeshInterface& h, const Vertex3d& v, double r) 
    : H(h), V(v), R(r) {}

  void BuildCellList(std::vector<int>&);
};

/*************************************************************/

class CylinderRefiner 
{
  const MeshInterface&  H;
  const Vertex3d&  V;
  double R;
  int D;

  bool QuadInCylinder(int) const;

 public:

  CylinderRefiner(const MeshInterface& h, const Vertex3d& v, 
		  double r, int d)
    : H(h), V(v), R(r), D(d) {}

  void BuildCellList(std::vector<int>&);
};

/*************************************************************/

class BallRefiner 
{
  const MeshInterface&  H;
  const Vertex3d&  V;
  double R;
  
  bool QuadInBall(int) const;

 public:
  
  BallRefiner(const MeshInterface& h, const Vertex3d& v, 
	      double r)
    : H(h), V(v), R(r) {}

  void BuildCellList(std::vector<int>&);
};
}

/*************************************************************/

#endif
