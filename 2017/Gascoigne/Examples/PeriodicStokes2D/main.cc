/**
*
* Copyright (C) 2009 by the Gascoigne 3D authors
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


#include  "localstokes.h"
#include  "localloop.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{

  ParamFile paramfile("periodic.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  LocalStokesProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("stokes", &LPD);
    
  FunctionalContainer FC;
  LocalDragFunctional LDF;
  LocalLiftFunctional LLF;
  FC.AddFunctional("Drag", &LDF);
  FC.AddFunctional("Lift", &LLF);

  LocalLoop loop;

  loop.BasicInit(&paramfile,&PC,&FC);
  loop.run("stokes");

  return 0;
}
