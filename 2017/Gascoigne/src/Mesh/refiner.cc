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


#include "refiner.h"

using namespace std;

/************************************************************/

namespace Gascoigne
{
bool PointRefiner::VertexInQuad(int i)
{
  Vertex2d V0 = H.vertex2d(H.vertex_of_cell(i,0));
  Vertex2d V1 = H.vertex2d(H.vertex_of_cell(i,1));
  Vertex2d V2 = H.vertex2d(H.vertex_of_cell(i,2));
  Vertex2d V3 = H.vertex2d(H.vertex_of_cell(i,3));

  A(0,0) = V1.x() - V0.x();
  A(0,1) = V3.x() - V0.x();
  A(1,0) = V1.y() - V0.y();
  A(1,1) = V3.y() - V0.y();

  A.gauss_jordan();

  DoubleVector y(2), z(2);
  y[0] = V.x()-V0.x();
  y[1] = V.y()-V0.y();

  A.multeq(z,y);

  int c=0;
  if (z[0]>=0) c++;
  if (z[0]<=1) c++;
  if (z[1]>=0) c++;
  if (z[1]<=1) c++;

  if (c==4) return 1;
  return 0;
}

/************************************************************/

void PointRefiner::BuildCellList(vector<int>& dst)
{
  for (int i=0; i<H.ncells(); i++)
    {
      if (VertexInQuad(i))
	{
	  dst.push_back(i);
	}
    }
}

/************************************************************/

bool CircleRefiner::QuadOnRadius(int q) const
{
  int out=0;
  int in=0;
  for (int i=0; i<8; i++)
    {
      Vertex3d V0 = H.vertex3d(H.vertex_of_cell(q,i));
      
      double dx = V.x()-V0.x();
      double dy = V.y()-V0.y();
      double dz = V.z()-V0.z();
      double d  = dx*dx+dy*dy+dz*dz - R*R;

      if (d<0) out++;
      else     in++;
    }
  if ((out==0) || (in==0)) return 0;
  
  return 1;
}

/************************************************************/

void CircleRefiner::BuildCellList(vector<int>& dst)
{
  for (int i=0; i<H.ncells(); i++)
    {
      if (QuadOnRadius(i))
	{
	  dst.push_back(i);
	}	  
    }
}


/************************************************************/

bool CylinderRefiner::QuadInCylinder(int q) const
{
  int i=0;
  double di[3];
  for (i=0; i<8; i++)
    {
      Vertex3d V0 = H.vertex3d(H.vertex_of_cell(q,i));
      
      di[0] = V.x()-V0.x();
      di[1] = V.y()-V0.y();
      di[2] = V.z()-V0.z();
      di[D]=0;
      double d  = di[0]*di[0]+di[1]*di[1]+di[2]*di[2] - R*R;
      
      if (d<0) break;
    }
  if (i==8) return 0;
  return 1;
}

/************************************************************/

void CylinderRefiner::BuildCellList(vector<int>& dst)
{
  for (int i=0; i<H.ncells(); i++)
    {
      if (QuadInCylinder(i))
	dst.push_back(i);
    }
  cerr << endl;
}

/************************************************************/

bool BallRefiner::QuadInBall(int q) const
{
  int i=0;
  double di[3];
  for (i=0; i<8; i++)
    {
      Vertex3d V0 = H.vertex3d(H.vertex_of_cell(q,i));
      
      di[0] = V.x()-V0.x();
      di[1] = V.y()-V0.y();
      di[2] = V.z()-V0.z();
      double d  = di[0]*di[0]+di[1]*di[1]+di[2]*di[2] - R*R;
      
      if (d<0) break;
    }
  if (i==8) return 0;
  return 1;
}

/************************************************************/

void BallRefiner::BuildCellList(vector<int>& dst)
{
  for (int i=0; i<H.ncells(); i++)
    {
      if (QuadInBall(i))
	dst.push_back(i);
    }
  cerr << endl;
}
}
