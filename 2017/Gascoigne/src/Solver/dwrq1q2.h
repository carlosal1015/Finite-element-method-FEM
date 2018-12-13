/**
*
* Copyright (C) 2004, 2005, 2008, 2009 by the Gascoigne 3D authors
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


#ifndef __DwrQ1Q2_h

#include "solverinterface.h"

/*-------------------------------------------------------*/

namespace Gascoigne
{
class DwrQ1Q2
{
 protected:

  SolverInterface& S;
  const ProblemDescriptorInterface* primalproblem;
  DiscretizationInterface*         discretization;

  virtual DiscretizationInterface* CreateOtherDiscretization() const;

  double ScalarProduct(nvector<double>& eta, const GlobalVector& f, 
		       const GlobalVector& z) const;

  double ScalarProduct(nvector<double>& eta, const VectorInterface& gf, 
		       const VectorInterface& gz) const;

 public:

  DwrQ1Q2(SolverInterface& SR);
  virtual ~DwrQ1Q2() {};

  double ScalarProductWithFluctuations(nvector<double>& eta, 
				       const VectorInterface& gf, 
				       const VectorInterface& gz) const;

  void PrimalResidualsHigher(VectorInterface& gf, 
			     const VectorInterface& gu);

  void DualResidualsHigher(VectorInterface& gf, const VectorInterface& gu, 
			   const VectorInterface& gz, 
			   const ProblemDescriptorInterface& PDI);

  double Estimator(nvector<double>& eta, VectorInterface& gf, 
		   const VectorInterface& gu, const VectorInterface& gz,
		   const ProblemDescriptorInterface& PDI);
  double EstimatorEnergy(nvector<double>& eta, VectorInterface& gf, 
			 const VectorInterface& gu);
};
}
/*-------------------------------------------------------*/

#endif
