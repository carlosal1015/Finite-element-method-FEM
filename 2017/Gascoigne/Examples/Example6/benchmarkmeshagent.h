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


#ifndef  __BenchmarkMeshAgent_h
#define  __BenchmarkMeshAgent_h

#include  "meshagent.h"

/* ----------------------------------------- */

class RunderKreis : public Gascoigne::BoundaryFunction<2>
{
  double              _r;
  Gascoigne::Vertex2d _c;

public :

  std::string GetName() const { return "RunderKreis";}
  void BasicInit(Gascoigne::Vertex2d c, double r) 
  {
      _c = c;
      _r = r;
    }
  double operator()(const Gascoigne::Vertex2d& c) const 
  {
      double r = - _r;
      for (int i=0; i<2; i++)
        {
          double dx = c[i]-_c[i];
          r += dx * dx;
        }
      return r;
    }
};

/*----------------------------------------------------------------------------*/

class CurvedMeshAgent : public Gascoigne::MeshAgent
{
 protected:

  RunderKreis RK;

 public:

  CurvedMeshAgent() : MeshAgent()
    {
      double r = 0.25;
      Gascoigne::Vertex2d v(2.,2.);
      RK.BasicInit(v,r);

      AddShape(80,&RK);
    }
};

/*-------------------------------------------------------------*/

#endif
