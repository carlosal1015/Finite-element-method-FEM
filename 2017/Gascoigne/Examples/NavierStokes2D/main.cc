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


#include  "local.h"


using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC) 
    {
      GetMeshAgentPointer() = new BenchMarkMeshAgent;
      StdLoop::BasicInit(paramfile, PC, FC);
    }
};

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("bench.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);
  
  ProblemContainer PC;
  PC.AddProblem("navier stokes", &LPD);

	    // Functionals
  FunctionalContainer FC;
  LocalDomainFunctionals_FlagForce df_flagforce_lift("lift");
  LocalDomainFunctionals_FlagForce df_flagforce_drag("drag");
  FC.AddFunctional("lift",&df_flagforce_lift);
  FC.AddFunctional("drag",&df_flagforce_drag);
  
  LocalLoop loop;
  loop.BasicInit(&paramfile, &PC, &FC);
  loop.run("navier stokes");

  return 0;
}
