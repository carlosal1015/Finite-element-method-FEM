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


#ifndef __EnergyEstimator_h

#include "solverinterface.h"
#include "q1.h"

/*-------------------------------------------------------*/

namespace Gascoigne
{

class EnergyEstimator
{
 protected:

  SolverInterface& S;
  const ProblemDescriptorInterface* primalproblem;
  Q1*                               discretization;
  double                            _d_visc;
  std::string                       _s_energytype;

 public:

  EnergyEstimator(SolverInterface& SR);
  ~EnergyEstimator() {};

  double Estimator(DoubleVector& eta, VectorInterface& gu, 
		   const VectorInterface& gf);
};

/*-------------------------------------------------------*/

}
#endif
