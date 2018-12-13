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


#ifndef  __LocalLoop_h
#define  __LocalLoop_h

#include  "stdloop.h"

#include  "meshagent.h"

/*-----------------------------------------*/

namespace Gascoigne
{

  class LocalLoop : public StdLoop
    {
    public:
      void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC=NULL)
      {
        GetMeshAgentPointer() = new MeshAgent;

        GetMeshAgent()->AddPeriodicMapping(90,92, new StdPeriodicMapping);
        GetMeshAgent()->AddPeriodicMapping(92,90, new StdPeriodicMapping);

        StdLoop::BasicInit(paramfile, PC, FC);
      };

    };
  
}

/*-----------------------------------------*/

#endif // __LocalLoop_h
