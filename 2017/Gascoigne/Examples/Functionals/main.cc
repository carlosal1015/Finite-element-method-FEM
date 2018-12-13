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


#include  "local1.h"
#include  "stdloop.h"

#include  "weightedpointfunctional.h"

using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh1.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("func", &LPD);
  
  /////////////
  // Functionals
  /////////////
    
  WeightedPointFunctional j0;
  vector<Vertex2d> v2d(2);
  v2d[0].x() = 0.25;
  v2d[0].y() = 0.25;
  
  v2d[1].x() = 0.5;
  v2d[1].y() = 0.5;    
  
  vector<int> comps(2,0);
  vector<double> w(2,1);

  j0.BasicInit(v2d,comps,w);
  j0.ExactValue()=0.09765625;

  LocalPointFunctional j1;
  j1.BasicInit(v2d,comps);

  LocalDragFunctional   j2; 
  LocalDomainFunctional j3;

  FunctionalContainer FC;
  FC.AddFunctional("0",&j0);
  FC.AddFunctional("1",&j1);
  FC.AddFunctional("2",&j2);
  FC.AddFunctional("3",&j3);
  
  
  /////////////
  // Loop
  /////////////
  StdLoop loop;
  loop.BasicInit(&paramfile, &PC, &FC);

  
  loop.run("func");

  return 0;
}
