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

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("box.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("laplace3d", &LPD);

  FunctionalContainer FC;
  LocalDragFunctional   j0; 
  LocalDomainFunctional j1;
  FC.AddFunctional("drag", &j0);
  FC.AddFunctional("domain", &j1);
  
  
  StdLoop loop;

  loop.BasicInit(&paramfile, &PC, &FC);

  loop.run("laplace3d");

  return 0;
}
