/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#ifndef  __StdTimeLoop_h
#define  __StdTimeLoop_h

#include  "stdloop.h"
#include  "timeinfo.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class StdTimeLoop : public virtual StdLoop
{
protected:

  TimeInfo    _timeinfo;
  virtual std::string SolveTimePrimal(VectorInterface& u, VectorInterface& f);

  virtual void TimeInfoBroadcast();
  void InitSolution(VectorInterface& u);

public:

  StdTimeLoop() : StdLoop() {}

  void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC,
		 const FunctionalContainer* FC=NULL);

  void run(const std::string& problemlabel);
  void adaptive_run(const std::string& problemlabel);
};
}

#endif
