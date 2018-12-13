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


#ifndef  __NonstationaryAlgorithm_h
#define  __NonstationaryAlgorithm_h

#include  "multilevelalgorithm.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class NonstationaryAlgorithm : public MultiLevelAlgorithm
{
 protected:

  double  dt, time, theta;

  void InitSolution(const std::string& initial, Gascoigne::VectorInterface& u) const;
  void TimeInfoBroadcast();

public:

  NonstationaryAlgorithm() :  MultiLevelAlgorithm() {}
  virtual ~NonstationaryAlgorithm() {}

  virtual void BasicInit(const Gascoigne::ParamFile* paramfile, MultiLevelSolver* MLS,
			 const Gascoigne::NumericInterface* NI,
			 const Gascoigne::ProblemContainer* PC);

  void ImplicitEuler            (const std::string&);
  void CrankNicholson           (const std::string&);
  void ThetaScheme              (const std::string&);
  void FractionalStepThetaScheme(const std::string&);
};
}

/*-----------------------------------------*/

#endif
